/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// ----------------------------------------------------------------------
// Contributing author: Ganga P Purja Pun (GMU)
// ----------------------------------------------------------------------


#include "pair_pinn.h"
#include <mpi.h>
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "utils.h"
#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

/*static const char cite_user_gmupinn_package[] =
  "USER-PINN package:\n\n"
  "@Article{Purja-Pun2020,\n"
  " author = {G. P. Purja Pun, V. Yamakov , J. Hickman , E. H. Glaessgen, Y. Mishin},\n"
  " title = {Development of a general-purpose machine-learning interatomic potential for aluminum by the physically informed neural network method },\n"
  " journal = {Phys. Rev. Materials},\n"
  " year =    2020,\n"
  " volume =  4,\n"
  " pages =   {113807}\n"
  "}\n\n";*/

#define MAXLINE 1024
#define DELTA 1.0e-20

// ----------------------------------------------------------------------

PairPINN::PairPINN(LAMMPS *lmp) : Pair(lmp)
{
  //if (lmp->citeme) lmp->citeme->add(cite_user_gmupinn_package);

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  // initialize arrays
  baseline_bo_params = nullptr;
  big_a = nullptr;
  alpha = nullptr;
  big_b = nullptr;
  beta = nullptr;
  small_h = nullptr;
  sigma = nullptr;
  small_a = nullptr;
  lambda = nullptr;

  deriv_wrt_big_a = nullptr;
  deriv_wrt_alpha = nullptr;
  deriv_wrt_big_b = nullptr;
  deriv_wrt_beta = nullptr;
  deriv_wrt_small_h = nullptr;
  deriv_wrt_sigma = nullptr;
  deriv_wrt_small_a = nullptr;
  deriv_wrt_lambda = nullptr;

}

// ----------------------------------------------------------------------

PairPINN::~PairPINN()
{
  if (copymode) return;

  if (params) {
    memory->destroy(params->mass);
    memory->destroy(params->LPOrders);
    memory->destroy(params->sigmas);

    if (params->layers) {
      for (int i=0; i<params->nlayers; i++) {
        memory->destroy (params->layers[i].Biases);
        memory->destroy (params->layers[i].Weights);
        memory->destroy (params->layers[i].fdot);
      }
      memory->destroy(params->layers);
    }
    delete params;
  }
  memory->destroy(elements);
  memory->destroy(map_gi);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
  }

  // free memory
  memory->destroy(baseline_bo_params);

  memory->destroy(big_a);
  memory->destroy(alpha);
  memory->destroy(big_b);
  memory->destroy(beta);
  memory->destroy(small_h);
  memory->destroy(sigma);
  memory->destroy(small_a);
  memory->destroy(lambda);

  memory->destroy(deriv_wrt_big_a);
  memory->destroy(deriv_wrt_alpha);
  memory->destroy(deriv_wrt_big_b);
  memory->destroy(deriv_wrt_beta);
  memory->destroy(deriv_wrt_small_h);
  memory->destroy(deriv_wrt_sigma);
  memory->destroy(deriv_wrt_small_a);
  memory->destroy(deriv_wrt_lambda);

  // delete accumulators
  memory->destroy(Sij);
  memory->destroy(Zij);
  memory->destroy(bij);
  memory->destroy(Sijk);
  memory->destroy(CSijk);
}

void PairPINN::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

void PairPINN::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  params = new Param();

  // read potential file and initialize potential parameters

  read_file(arg[2]);

  setup_params ();

  // create a map to j-k pairs
  // for a ternary system; the table should look like this:
  //
  // | 0  1  2 |
  // | 1  3  4 |
  // | 2  4  5 |

  map_gi = memory->create(map_gi, nelements, nelements, "create:map_gi");

  int tmp = 0;
  for (int i=0; i<nelements; i++) {
    for (int j=i; j<nelements; j++) {
      map_gi[i][j] = map_gi[j][i] = tmp;
      tmp++;
    }
  }


  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i], "NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i], elements[j]) == 0) break;
    if (j < nelements) map[i-2] = j;
    else error->all(FLERR,"No matching element in PINN potential file");
  }

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR,i,params->mass[map[i]]);
        count++;
      }
      //scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

void PairPINN::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style PINN requires newton pair on");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

// ----------------------------------------------------------------------

void PairPINN::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  map = new int[n+1]; // map for atom types
}

// ----------------------------------------------------------------------

void PairPINN::setup_params()
{
  // allocate memory for local BO parameter related arrays
  // sizes of the arrays do not change once they are allocated here

  if (nelements <= 0) {
    error->all (FLERR, fmt::format("Illegal # of unique elements: {}",
                                   nelements));
  }

  big_a = memory->create(big_a, nelements, nelements, "BOP set: A[i,j]");
  alpha = memory->create(alpha, nelements, nelements, "BOP set: alpha[i,j]");
  big_b = memory->create(big_b, nelements, nelements, "BOP set: B[i,j]");
  beta = memory->create(beta, nelements, nelements, "BOP set: beta[i,j]");
  small_h = memory->create(small_h, nelements, nelements, nelements,
                           "BOP set: h[i,j,k]");
  sigma = memory->create(sigma, nelements, "BOP set: sigma[i]");
  small_a = memory->create(small_a, nelements, nelements, nelements,
                           "BOP set: a[i,j,k]");
  lambda = memory->create(lambda, nelements, nelements, nelements,
                          "BOP set: lambda[i,j,k]");

  deriv_wrt_big_a = memory->create(deriv_wrt_big_a, nelements, nelements,
                                   "create:init_one");
  deriv_wrt_alpha = memory->create(deriv_wrt_alpha, nelements, nelements,
                                   "create:init_one");
  deriv_wrt_big_b = memory->create(deriv_wrt_big_b, nelements, nelements,
                                   "create:init_one");
  deriv_wrt_beta = memory->create(deriv_wrt_beta, nelements, nelements,
                                  "create:init_one");
  deriv_wrt_small_h = memory->create(deriv_wrt_small_h, nelements, nelements,
                                     nelements, "create:init_one");
  deriv_wrt_sigma = memory->create(deriv_wrt_sigma, nelements,
                                   "create:init_one");
  deriv_wrt_small_a = memory->create(deriv_wrt_small_a, nelements, nelements,
                                     nelements, "create:init_one");
  deriv_wrt_lambda = memory->create(deriv_wrt_lambda, nelements, nelements,
                                    nelements, "create:init_one");

  // set cutmax to max of all params
  cutmax = 1.5 * params->cut;
  cutmaxsq = cutmax * cutmax;
}


// ----------------------------------------------------------------------

double PairPINN::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

// ----------------------------------------------------------------------

void PairPINN::compute(int eflag, int vflag)
{
  int nblock, unit_block, ngi;
  int i, j, k, l, ii, jj, kk, ll, inum, jnum, knum;
  int itype, jtype, ktype, ltype;
  double fpair, evdwl;
  tagint itag,jtag;
  int *ilist, *jlist, *klist, *numneigh, **firstneigh;
  double *gilist;
  Layer *atomnnlayers;
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double xij, yij, zij, drij[3];
  double xik, yik, zik, drik[3];
  double xjk, yjk, zjk, drjk[3];
  double rsqij, rsqik, rsqjk;
  double r, rij, rik, rjk;
  double xtmp, ytmp, ztmp;
  double *Pl, *Pl_d;
  double Ep;
  double dbij;

  // initialize the accumulators

  Sij = nullptr;
  Zij = nullptr;
  bij = nullptr;
  Sijk = nullptr;
  CSijk = nullptr;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // size of feature vecotr
  nblock = nelements * (nelements + 1) / 2;
  unit_block = params->nsigmas * params->nLPOrders;
  ngi = pinn_kind_flag1 * nblock * unit_block;

  // allocate memory for Legendre polynomial orders
  Pl = new double [params->LPOrders[params->nLPOrders - 1] + 1];
  Pl[0] = 1.0; // order l = 0

  // allocate memory for derivatives of Legendre polynomial orders
  Pl_d = new double [params->LPOrders[params->nLPOrders - 1] + 1];
  Pl_d[0] = 0.0; // derivative of order l = 0
  Pl_d[1] = 1.0; // derivative of order l = 1

  gilist = memory->create(gilist, ngi, "compute:gilist");

  atomnnlayers = memory->create(atomnnlayers, params->nlayers,
                                "compute:atomnnlayers");

  for (j = 1; j < params->nlayers; j++) { // loop over NN layers
    atomnnlayers[j].fdot = nullptr;
    atomnnlayers[j].fdot = memory->create(atomnnlayers[j].fdot,
                                          params->layers[j].nnodes,
                                          "compute:atomnnlayers");
  }

  // allocate memory for feature vector

  double *Gx = new double[ngi];
  double *Gy = new double[ngi];
  double *Gz = new double[ngi];

  // allocate memory for actual baseline BOPs and perturbations

  double *actual_bops = new double[nbaseline_bo_params];
  double *perturb_bops = new double[nbaseline_bo_params];
  double *array_derivs_wrt_bops = new double[nbaseline_bo_params];

  // X, Y and Z components of derivatives of local BOPs

  double *nnfx = new double[nbaseline_bo_params];
  double *nnfy = new double[nbaseline_bo_params];
  double *nnfz = new double[nbaseline_bo_params];

  // compute forces due to Gis through NN

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = map[type[i]];

    // initialize feature vector
    for (int q = 0; q < ngi; q++)
      gilist[q] = params->giref;

    // compute feature vector
    compute_atomic_LSP(i, gilist); // transformed Gis !!!

    // initialize actual BOPs
    memset(actual_bops, 0, nbaseline_bo_params * sizeof (double));

    // run through NN
    eval_nnet(gilist, ngi, perturb_bops,
              params->layers[params->nlayers - 1].nnodes);

    // compute actual BOPs
    for (int q = 0; q < nbaseline_bo_params; ++q) {
      if (nBase) actual_bops[q] = baseline_bo_params[q] + perturb_bops[q];
      else actual_bops[q] = perturb_bops[q];
    }

    // partition elements of array of BOPs into respective BO parameter sets
    unpack_vec_to_bop_sets (actual_bops, nbaseline_bo_params);

    // copy fdot from NN layers to atomnnlayers for later use
    for (j = 1; j < params->nlayers-1; j++) { // safe to use 'j'
      for (k = 0; k < params->layers[j].nnodes; k++) // safe to use 'k'
        atomnnlayers[j].fdot[k] = params->layers[j].fdot[k];
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jnum = numneigh[i];
    jlist = firstneigh[i];
    knum = numneigh[i];
    klist = firstneigh[i];

    // ------------------------------------------------------------------
    //
    // compute E_i and its derivatives wrt 'r_i' at constant BOPs
    //
    // ------------------------------------------------------------------

    Sij = memory->create(Sij, jnum, "grow:Sij");
    Zij = memory->create(Zij, jnum, "grow:Zij");
    bij = memory->create(bij, jnum, "grow:bij");
    Sijk = memory->create(Sijk, jnum, jnum, "grow:Sijk");
    CSijk = memory->create(CSijk, jnum, jnum, "grow:CSijk");

    // initialize Sij, Zij, bij, Sijk and CSijk

    for (int m = 0; m < jnum; m++) {
      Sij[m] = 1.0;
      Zij[m] = 0.0;
      bij[m] = 0.0;
      for (int k = 0; k < jnum; k++) {
        Sijk[m][k] = 1.0;
        CSijk[m][k] = 1.0;
      }
    }

    // two-body interactions

    for (jj = 0; jj < jnum; jj++) { // loop over all neighbors of i
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];

      drij[0] = x[j][0] - xtmp;
      drij[1] = x[j][1] - ytmp;
      drij[2] = x[j][2] - ztmp;

      rsqij = dot_product (drij, drij, 3);

      if (rsqij > params->cutsq) continue;

      rij = sqrt(rsqij);

      // repulsive term

      if (eflag) evdwl = 0.5 * bop_fr(itype, jtype, rij) * ann_fc(rij, params);

      fpair = 0.5 * (ann_fc_d (rij, params) - ann_fc (rij, params) *
                     alpha[itype][jtype]) * bop_fr (itype, jtype, rij) / rij;

      f[i][0] += drij[0] * fpair;    // convention is r_ij = r_j - r_i
      f[i][1] += drij[1] * fpair;
      f[i][2] += drij[2] * fpair;
      f[j][0] -= drij[0] * fpair;
      f[j][1] -= drij[1] * fpair;
      f[j][2] -= drij[2] * fpair;

      if (evflag) {
        ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, drij[0],
            drij[1], drij[2]);
      }

      // compute S_ijk and S_ij

      for (kk = 0; kk < knum; kk++) { // loop over all neighbors of i
        if (jj == kk) continue;
        k = klist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];

        drik[0] = x[k][0] - xtmp;
        drik[1] = x[k][1] - ytmp;
        drik[2] = x[k][2] - ztmp;

        rsqik = dot_product (drik, drik, 3);
        rik = sqrt(rsqik);

        drjk[0] = drik[0] - drij[0];
        drjk[1] = drik[1] - drij[1];
        drjk[2] = drik[2] - drij[2];

        rsqjk = dot_product (drjk, drjk, 3);
        rjk = sqrt(rsqjk);
        r = rik + rjk - rij;

        if (r < params->cut) {
          Sijk[jj][kk] = 1.0 - ann_fc (r, params) *
              exp(-lambda[itype][jtype][ktype] * fabs(r));
          Sij[jj] *= Sijk[jj][kk];
        }
      }
    }

    // three-body interactions

    Ep = 0.0;

    for (jj = 0; jj < jnum; jj++) { // loop over all neighbors of i
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];

      drij[0] = x[j][0] - xtmp;
      drij[1] = x[j][1] - ytmp;
      drij[2] = x[j][2] - ztmp;

      rsqij = dot_product (drij, drij, 3);

      if (rsqij > params->cutsq) continue;

      rij = sqrt(rsqij);

      // compute cosine of theta_ijk and z_ij

      for (kk = 0; kk < knum; kk++) {
        if (jj == kk) continue;
        k = klist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];

        drik[0] = x[k][0] - xtmp;
        drik[1] = x[k][1] - ytmp;
        drik[2] = x[k][2] - ztmp;

        rsqik = dot_product (drik, drik, 3);

        if (rsqik > params->cutsq) continue;

        rik = sqrt(rsqik);

        CSijk[jj][kk] = dot_product (drij, drik, 3) / (rij * rik);
        Zij[jj] += bop_zeta (itype, jtype, ktype, params, Sij[kk],
                             CSijk[jj][kk], rik);
      }

      if (Zij[jj] <= -1.0) {
        printf("rank:%d atom:%d --> Zij[%d] = %e\n", comm->me, i, j, Zij[jj]);
        fflush(stdout);
        bij[jj] = 1.0 / sqrt(fabs(1.0 + Zij[jj]) + DELTA);
      } else {
        bij[jj] = 1.0 / sqrt(1.0 + Zij[jj]);
      }

      dbij = cube(bij[jj]) * copysign (1.0, 1.0 + Zij[jj]);

      // accumulate promotional energy

      Ep += bij[jj] * Sij[jj] * ann_fc (rij, params);

      // attractive term

      if (eflag) {
        evdwl = -0.5 * Sij[jj] * bij[jj] * ann_fc(rij, params) *
            exp(big_b[itype][jtype] - beta[itype][jtype] * rij);
      }

      fpair = 0.5 * (beta[itype][jtype] * ann_fc(rij, params)
                     - ann_fc_d(rij, params)) * Sij[jj] * bij[jj] *
          exp(big_b[itype][jtype] - beta[itype][jtype] * rij) / rij;

      f[i][0] += drij[0] * fpair;    // convention is r_ij = r_j - r_i
      f[i][1] += drij[1] * fpair;
      f[i][2] += drij[2] * fpair;
      f[j][0] -= drij[0] * fpair;
      f[j][1] -= drij[1] * fpair;
      f[j][2] -= drij[2] * fpair;

      if (evflag) {
        ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, drij[0],
            drij[1], drij[2]);
      }

      for (kk = 0; kk < knum; kk++) { // loop over all neighbors of i
        if (jj == kk) continue;
        k = klist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];

        drik[0] = x[k][0] - xtmp;
        drik[1] = x[k][1] - ytmp;
        drik[2] = x[k][2] - ztmp;

        rsqik = dot_product (drik, drik, 3);

        rik = sqrt(rsqik);

        if (rsqik < params->cutsq) {
          // derivative of b_ij wrt 'r_i' from attractive term
          fpair = 0.5 * Sij[jj] * ann_fc(rij, params)
              * exp(big_b[itype][jtype] - beta[itype][jtype] * rij);

          // prefactor for derivative of cosine_theta_ijk wrt 'r_i'
          double prefactor_cs_ijk = fpair * -dbij *
              small_a[itype][jtype][ktype] * Sij[kk] * ann_fc(rik, params) *
              (CSijk[jj][kk] - small_h[itype][jtype][ktype]) ;

          // prefactor for derivative of f_c(r_ik) wrt 'r_i'
          fpair *= 0.5 * dbij * small_a[itype][jtype][ktype] * Sij[kk] *
              square(CSijk[jj][kk] - small_h[itype][jtype][ktype]) *
              ann_fc_d(rik, params) / rik;

          f[i][0] += drik[0] * fpair;    // convention is r_ij = r_j - r_i
          f[i][1] += drik[1] * fpair;
          f[i][2] += drik[2] * fpair;
          f[k][0] -= drik[0] * fpair;
          f[k][1] -= drik[1] * fpair;
          f[k][2] -= drik[2] * fpair;

          if (vflag_atom) {
            v_tally2(i, k, fpair, drik);
          }

          double fi[3], fj[3], fk[3];

          force_costheta_ijk_ri (prefactor_cs_ijk, CSijk[jj][kk], drij,
                                 drik, fi, fj, fk);

          f[i][0] += fi[0];
          f[i][1] += fi[1];
          f[i][2] += fi[2];
          f[j][0] -= fj[0];
          f[j][1] -= fj[1];
          f[j][2] -= fj[2];
          f[k][0] -= fk[0];
          f[k][1] -= fk[1];
          f[k][2] -= fk[2];

          if (vflag_atom) v_tally3(i, j, k, fj, fk, drij, drik);

          // derivative of S_ik wrt 'r_i'

          double dril[3], ril, drkl[3], rkl;

          for (ll = 0; ll < knum; ll++) {
            if (ll == kk) continue;
            l = klist[ll];
            l &= NEIGHMASK;
            ltype = map[type[l]];

            dril[0] = x[l][0] - xtmp;
            dril[1] = x[l][1] - ytmp;
            dril[2] = x[l][2] - ztmp;

            ril = sqrt(dot_product (dril, dril, 3));

            drkl[0] = dril[0] - drik[0];
            drkl[1] = dril[1] - drik[1];
            drkl[2] = dril[2] - drik[2];

            rkl = sqrt(dot_product (drkl, drkl, 3));

            r = ril + rkl - rik;

            if (r < params->cut) {
              double prefactor = -0.5 * Sij[jj] * bop_fa(itype, jtype, rij)
                  * ann_fc(rij, params);

              prefactor *= -0.5 * dbij * small_a[itype][jtype][ktype] *
                  square(CSijk[jj][kk] - small_h[itype][jtype][ktype])
                  * ann_fc(rik, params);

              prefactor *= Sij[kk] / Sijk[kk][ll] *
                  exp(-lambda[itype][ktype][ltype] * fabs(r)) *
                  (lambda[itype][ktype][ltype] * ann_fc (r, params)
                   - ann_fc_d(r, params));

              double fi[3], fk[3], fl[3];

              force_Sijk_ri (prefactor, drik, dril, rik, ril, rkl,
                             CSijk[kk][ll], fi, fk, fl);

              f[i][0] += fi[0];
              f[i][1] += fi[1];
              f[i][2] += fi[2];
              f[k][0] -= fk[0];
              f[k][1] -= fk[1];
              f[k][2] -= fk[2];
              f[l][0] -= fl[0];
              f[l][1] -= fl[1];
              f[l][2] -= fl[2];

              if (vflag_atom) v_tally3(i,k,l,fk,fl,drik,dril);
            }
          }
        }

        drjk[0] = drik[0] - drij[0];
        drjk[1] = drik[1] - drij[1];
        drjk[2] = drik[2] - drij[2];

        rsqjk = dot_product (drjk, drjk, 3);

        rjk = sqrt(rsqjk);

        r = rik + rjk - rij;

        if (r < params->cut) {
          double prefac = -0.5 * bij[jj] * bop_fa (itype, jtype, rij)
              * ann_fc(rij, params);

          prefac *= Sij[jj] / Sijk[jj][kk] *
              (lambda[itype][jtype][ktype] * ann_fc (r, params)
               - ann_fc_d (r, params))
              * exp(-lambda[itype][jtype][ktype] * fabs(r));

          double fi[3], fj[3], fk[3];

          force_Sijk_ri (prefac, drij, drik, rij, rik, rjk,
                         CSijk[jj][kk], fi, fj, fk);

          f[i][0] += fi[0];
          f[i][1] += fi[1];
          f[i][2] += fi[2];
          f[j][0] -= fj[0];
          f[j][1] -= fj[1];
          f[j][2] -= fj[2];
          f[k][0] -= fk[0];
          f[k][1] -= fk[1];
          f[k][2] -= fk[2];

          if (vflag_atom) v_tally3(i,j,k,fj,fk,drij,drik);
        }
      }
    }

    // update promotional energy

    if (evflag) ev_tally_full(i, -2.0 * sigma[itype] * sqrt(Ep),
                              0.0, 0.0, 0.0, 0.0, 0.0);

    // -------------------------------------------------------------
    //
    // forces from promotional energy
    //
    // -------------------------------------------------------------

    for (jj = 0; jj < jnum; jj++) { // loop over all neighbors of i
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];

      drij[0] = x[j][0] - xtmp;
      drij[1] = x[j][1] - ytmp;
      drij[2] = x[j][2] - ztmp;

      rsqij = dot_product (drij, drij, 3);

      if (rsqij > params->cutsq) continue;

      rij = sqrt(rsqij);

      dbij = cube(bij[jj]) * copysign (1.0, 1.0 + Zij[jj]);

      // derivative of f_c(r_ij) w.r.t. 'r_i'

      fpair = -0.5 * sigma[itype] / sqrt(Ep) * Sij[jj] * bij[jj]
          * ann_fc_d (rij, params) / rij;

      f[i][0] += drij[0] * fpair;    // convention is r_ij = r_j - r_i
      f[i][1] += drij[1] * fpair;
      f[i][2] += drij[2] * fpair;
      f[j][0] -= drij[0] * fpair;
      f[j][1] -= drij[1] * fpair;
      f[j][2] -= drij[2] * fpair;

      if (vflag_atom) v_tally2 (i, j, fpair, drij);

      // derivative of S_ij w.r.t. 'r_i'

      for (kk = 0; kk < knum; kk++) { // again loop over all neighbors of i
        if (jj == kk) continue;
        k = klist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];

        drik[0] = x[k][0] - xtmp;
        drik[1] = x[k][1] - ytmp;
        drik[2] = x[k][2] - ztmp;

        rsqik = dot_product (drik, drik, 3);
        rik = sqrt(rsqik);

        if (rsqik < params->cutsq) {
          // prefactor of derivative of b_ij w.r.t. 'r_i'
          fpair = 0.5 * sigma[itype] / sqrt(Ep) * Sij[jj]
              * ann_fc(rij, params);

          // prefactor of derivative of costheta_ijk w.r.t. 'r_i'
          double prefactor_cs_ijk = fpair * -dbij
              * small_a[itype][jtype][ktype] * Sij[kk]
              * (CSijk[jj][kk] - small_h[itype][jtype][ktype])
              * ann_fc (rik, params);

          // prefactor of derivative of f_c(r_ik) w.r.t. 'r_i'
          fpair *= 0.5 * dbij * small_a[itype][jtype][ktype] * Sij[kk]
              * square(CSijk[jj][kk] - small_h[itype][jtype][ktype])
              * ann_fc_d(rik, params) / rik;

          f[i][0] += drik[0] * fpair;    // convention is r_ij = r_j - r_i
          f[i][1] += drik[1] * fpair;
          f[i][2] += drik[2] * fpair;
          f[k][0] -= drik[0] * fpair;
          f[k][1] -= drik[1] * fpair;
          f[k][2] -= drik[2] * fpair;


          if (vflag_atom) v_tally2 (i, k, fpair, drik);

          // derivative of cosine_theta w.r.t. 'r_i'

          double fi[3], fj[3], fk[3];

          force_costheta_ijk_ri (prefactor_cs_ijk, CSijk[jj][kk], drij, drik,
                                 fi, fj, fk);

          f[i][0] += fi[0];
          f[i][1] += fi[1];
          f[i][2] += fi[2];
          f[j][0] -= fj[0];
          f[j][1] -= fj[1];
          f[j][2] -= fj[2];
          f[k][0] -= fk[0];
          f[k][1] -= fk[1];
          f[k][2] -= fk[2];

          if (vflag_atom) v_tally3(i,j,k,fj,fk,drij,drik);

          // derivative of S_ik wrt 'r_i'

          double dril[3], ril, drkl[3], rkl;

          for (ll = 0; ll < knum; ll++) {
            if (ll == kk) continue;
            l = klist[ll];
            l &= NEIGHMASK;
            ltype = map[type[l]];

            dril[0] = x[l][0] - xtmp;
            dril[1] = x[l][1] - ytmp;
            dril[2] = x[l][2] - ztmp;

            ril = sqrt(dot_product (dril, dril, 3));

            drkl[0] = dril[0] - drik[0];
            drkl[1] = dril[1] - drik[1];
            drkl[2] = dril[2] - drik[2];

            rkl = sqrt(dot_product (drkl, drkl, 3));

            r = ril + rkl - rik;

            if (r < params->cut) {
              // prefactor of derivative of S_ik w.r.t. 'r_i'
              double prefactor = 0.5 * sigma[itype] / sqrt(Ep) * Sij[jj]
                  * ann_fc (rij, params);

              prefactor *= -0.5 * dbij * small_a[itype][jtype][ktype] *
                  square(CSijk[jj][kk] - small_h[itype][jtype][ktype])
                  * ann_fc (rik, params);

              prefactor *= Sij[kk] / Sijk[kk][ll]
                  * exp(-lambda[itype][ktype][ltype] * r)
                  * (lambda[itype][ktype][ltype] * ann_fc (r, params)
                     - ann_fc_d (r, params));

              double fi[3], fk[3], fl[3];

              force_Sijk_ri (prefactor, drik, dril, rik, ril, rkl,
                             CSijk[kk][ll], fi, fk, fl);

              f[i][0] += fi[0];
              f[i][1] += fi[1];
              f[i][2] += fi[2];
              f[k][0] -= fk[0];
              f[k][1] -= fk[1];
              f[k][2] -= fk[2];
              f[l][0] -= fl[0];
              f[l][1] -= fl[1];
              f[l][2] -= fl[2];

              if (vflag_atom) v_tally3(i, k, l, fk, fl, drik, dril);
            }
          }
        }

        drjk[0] = drik[0] - drij[0];
        drjk[1] = drik[1] - drij[1];
        drjk[2] = drik[2] - drij[2];

        rsqjk = dot_product (drjk, drjk, 3);

        rjk = sqrt(rsqjk);

        r = rik + rjk - rij;

        if (r < params->cut) {
          // prefactor of derivative of S_ij w.r.t. 'r_i'
          double prefac = 0.5 * sigma[itype] / sqrt(Ep) * bij[jj]
              * ann_fc (rij, params);

          prefac *= Sij[jj] / Sijk[jj][kk]
              * exp(-lambda[itype][jtype][ktype] * r)
              * (lambda[itype][jtype][ktype] * ann_fc (r, params)
                 - ann_fc_d (r, params));

          double fi[3], fj[3], fk[3];

          force_Sijk_ri (prefac, drij, drik, rij, rik, rjk, CSijk[jj][kk],
                         fi, fj, fk);

          f[i][0] += fi[0];
          f[i][1] += fi[1];
          f[i][2] += fi[2];
          f[j][0] -= fj[0];
          f[j][1] -= fj[1];
          f[j][2] -= fj[2];
          f[k][0] -= fk[0];
          f[k][1] -= fk[1];
          f[k][2] -= fk[2];

          if (vflag_atom) v_tally3(i, j, k, fj, fk, drij, drik);
        }
      }
    }

    // --------------------------------------------------------------
    //
    // compute derivatives of E_i wrt to BOPs
    //
    iderivs_wrt_bops (i, sqrt(Ep));
    //
    // pack derivatives wrt to BOPs
    pack_derivs_wrt_bop_sets_to_vec (array_derivs_wrt_bops,
                                     nbaseline_bo_params);
    //
    // --------------------------------------------------------------


    // --------------------------------------------------------------
    //
    // compute derivatives of local BOPs wrt 'r_i'
    //

    for (jj = 0; jj < jnum; jj++) { // loop over all neighbors of i
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];

      drij[0] = x[j][0] - xtmp;
      drij[1] = x[j][1] - ytmp;
      drij[2] = x[j][2] - ztmp;

      rsqij = dot_product (drij, drij, 3);

      // initialize feature vector
      for (int n = 0; n < ngi; ++n) {
        Gx[n] = 0.0;
        Gy[n] = 0.0;
        Gz[n] = 0.0;
      }

      // if neighbor j is within i's cutoff radius

      if (rsqij > cutmaxsq) continue;

      rij = sqrt(rsqij);

      for (kk = 0; kk < knum; kk++) { // loop over all neighbors of i again
        k = klist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];

        drik[0] = x[k][0] - xtmp;
        drik[1] = x[k][1] - ytmp;
        drik[2] = x[k][2] - ztmp;
        rsqik = dot_product (drik, drik, 3);

        if (rsqik > cutmaxsq) continue;
        rik = sqrt(rsqik);

        double cosijk = dot_product (drij, drik, 3) / (rik * rij);
        Pl[1] = cosijk;
        double dcosijk_X = (drik[0] / rik - drij[0] / rij * cosijk) / rij;
        double dcosijk_Y = (drik[1] / rik - drij[1] / rij * cosijk) / rij;
        double dcosijk_Z = (drik[2] / rik - drij[2] / rij * cosijk) / rij;

        // compute Legendre polynomials and their
        // corresponding derivatives over ijk
        for (int l=1; l<params->LPOrders[params->nLPOrders - 1]; l++) {
          Pl[l + 1] = ((2.0 * l + 1.0) * cosijk * Pl[l] - l
                       * Pl[l - 1]) / (l + 1);
          Pl_d[l + 1] = (l + 1) * Pl[l] + cosijk * Pl_d[l];
        }

        for(int m = 0; m<params->nsigmas; m++) { // loop over r0s
          double fs_rij = pinn_fs (rij, m);
          double fs_rik  = pinn_fs (rik, m);
          double fs_rij_d = pinn_fs_d (rij, m);

          // loop over Legendre polynomial orders
          for(int l = 0; l<params->nLPOrders; l++) {
            double dGx = Pl[params->LPOrders[l]] * fs_rij_d * (drij[0] / rij)
                * fs_rik + Pl_d[params->LPOrders[l]] * dcosijk_X * fs_rij
                * fs_rik;
            double dGy = Pl[params->LPOrders[l]] * fs_rij_d * (drij[1] / rij)
                * fs_rik + Pl_d[params->LPOrders[l]] * dcosijk_Y * fs_rij
                * fs_rik;
            double dGz = Pl[params->LPOrders[l]] * fs_rij_d * (drij[2] / rij)
                * fs_rik + Pl_d[params->LPOrders[l]] * dcosijk_Z * fs_rij
                * fs_rik;

            // update G-vector
            int gindx = itype * nblock * unit_block * pinn_kind_flag2
                + map_gi[jtype][ktype] * unit_block + l * params->nsigmas + m;

            // use the original Gis not the transformed here !!!
            double denom = sqrt(square(sinh(gilist[gindx])) + 1.0);

            Gx[gindx] += 2.0 * dGx / denom;
            Gy[gindx] += 2.0 * dGy / denom;
            Gz[gindx] += 2.0 * dGz / denom;
          }
          // -------------- end of loop over Legendre polynomials -----------
        }
        // ------------- end of loop over r0s -------------
      }
      // ------------ end of loop over kk neighbors of i -----------

      // run through ANN for i-j bond
      eval_nnet_d(atomnnlayers, Gx, params->layers[0].nnodes, nnfx,
          params->layers[params->nlayers - 1].nnodes);
      eval_nnet_d(atomnnlayers, Gy, params->layers[0].nnodes, nnfy,
          params->layers[params->nlayers - 1].nnodes);
      eval_nnet_d(atomnnlayers, Gz, params->layers[0].nnodes, nnfz,
          params->layers[params->nlayers - 1].nnodes);

      double fx = dot_product (array_derivs_wrt_bops, nnfx,
                               nbaseline_bo_params);
      double fy = dot_product (array_derivs_wrt_bops, nnfy,
                               nbaseline_bo_params);
      double fz = dot_product (array_derivs_wrt_bops, nnfz,
                               nbaseline_bo_params);

      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;

      f[j][0] -= fx;
      f[j][1] -= fy;
      f[j][2] -= fz;

      if (vflag_atom) {
        ev_tally_xyz (i, j, nlocal, newton_pair, 0.0, 0.0, fx, fy, fz,
                      drij[0], drij[1], drij[2]);
      }
      // ---- end of force computation on atom i due to neighbors rij ----
    }
    // end of loop over jj neighbors of i
    //
    // --------------------------------------------------------------

    memory->destroy(Sij);
    memory->destroy(Zij);
    memory->destroy(bij);
    memory->destroy(Sijk);
    memory->destroy(CSijk);
  }
  // ----------- end of loop over my atoms i ---------

  //clear memory
  memory->destroy(gilist);
  memory->destroy(atomnnlayers);

  if (Pl) delete [] Pl;
  if (Pl_d) delete [] Pl_d;
  if (Gx) delete [] Gx;
  if (Gy) delete [] Gy;
  if (Gz) delete [] Gz;

  if (actual_bops) delete [] actual_bops;
  if (perturb_bops) delete [] perturb_bops;
  if (array_derivs_wrt_bops) delete [] array_derivs_wrt_bops;
  if (nnfx) delete [] nnfx;
  if (nnfy) delete [] nnfy;
  if (nnfz) delete [] nnfz;

  // compute pressure via f dot r
  if (vflag_fdotr) virial_fdotr_compute();
}

void PairPINN::pack_derivs_wrt_bop_sets_to_vec(double *pvec, const int n)
{
  // ------------------------------------------
  // Sets of derivatives wrt BOPs are packed
  // into a given array of size 'n'.
  // The order is important !!!
  //
  // Current order:
  //
  // A, alpha, B, beta, h, sigma, a, lambda
  // ------------------------------------------

  // check size of input vector
  if (n <= 0) error->all (FLERR, "Invalid size (n)");

  int pcount = 0;

  for (int n = 0; n < nelements; n++) { // loop over number of types
    // As
    for (int i = 0; i < nelements; i++) {
      pvec[pcount++] = deriv_wrt_big_a[n][i];
    }

    // alphas
    for (int i = 0; i < nelements; i++) {
      pvec[pcount++] = deriv_wrt_alpha[n][i];
    }

    // Bs
    for (int i = 0; i < nelements; i++) {
      pvec[pcount++] = deriv_wrt_big_b[n][i];
    }

    // betas
    for (int i = 0; i < nelements; i++) {
      pvec[pcount++] = deriv_wrt_beta[n][i];
    }

    // hs
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        pvec[pcount++] = deriv_wrt_small_h[n][i][j];
      }
    }

    // sigmas
    pvec[pcount++] = deriv_wrt_sigma[n];

    // a's
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        pvec[pcount++] = deriv_wrt_small_a[n][i][j];
      }
    }

    // lambdas
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        pvec[pcount++] = deriv_wrt_lambda[n][i][j];
      }
    }
  }

  if (pcount != n) {
    error->all (FLERR, "Insufficient number of elements read");
  }
}

// ----------------------------------------------------------------------

void PairPINN::force_Sijk_ri(double pref, double *dri1, double *dri2,
                             double ri1, double ri2, double r12,
                             double cstheta, double *fi, double *f1,
                             double *f2)
{
  dcostheta_ijk_dri (cstheta, dri1, dri2, f1, f2);

  for (int q = 0; q < 3; q++) {
    f1[q] = dri1[q] / ri1 - (dri1[q] - ri2 * cstheta * dri1[q] / ri1) / r12
        - ri1 * ri2 * f1[q] / r12;
    f1[q] *= pref;
    f2[q] = -dri2[q] / ri2 - (dri2[q] - ri1 * cstheta * dri2[q] / ri2) / r12
        - ri1 * ri2 * f2[q] / r12;
    f2[q] *= pref;
    fi[q] = f1[q] + f2[q];
  }
}

void PairPINN::dcostheta_ijk_dri(double cstheta, double *drj,
                                 double *drk, double *fj, double *fk)
{
  double rij, rik;

  rij = sqrt(dot_product (drj, drj, 3));
  rik = sqrt(dot_product (drk, drk, 3));

  for (int i = 0; i < 3; i++) {
    fj[i] = (-drk[i] / rik + cstheta * drj[i] / rij) / rij;
    fk[i] = (-drj[i] / rij + cstheta * drk[i] / rik) / rik;
  }
}

void PairPINN::force_costheta_ijk_ri(double pref, double cstheta,
                                     double *drj, double *drk,
                                     double *fi, double *fj, double *fk)
{
  double rij, rik;

  rij = sqrt(dot_product (drj, drj, 3));
  rik = sqrt(dot_product (drk, drk, 3));

  for (int i = 0; i < 3; i++) {
    fj[i] = pref * (-drk[i] / rik + cstheta * drj[i] / rij) / rij;
    fk[i] = pref * (-drj[i] / rij + cstheta * drk[i] / rik) / rik;
    fi[i] = fj[i] + fk[i];
  }
}

double PairPINN::dot_product(double *x, double *y, int n)
{
  double tmp = 0.0;

  for (int i = 0; i < n; ++i) {
    tmp += x[i] * y[i];
  }

  return tmp;
}


// ----------------------------------------------------------------------

double PairPINN::bop_fa(int it, int jt, double r)
{
  // only -f_A(r)
  return -exp(big_b[it][jt] - beta[it][jt] * r);
}

double PairPINN::bop_fa_d(int it, int jt, Param *param, double r)
{
  // derivative of f_A(r)f_c(r) wrt 'r'

  double tmp_fc, tmp_fc_d, tmp_exp;

  tmp_fc = ann_fc (r, param);
  tmp_fc_d = ann_fc_d (r, param);
  tmp_exp = bop_fa (it, jt, r);

  return tmp_exp * (tmp_fc_d - tmp_fc * beta[it][jt]) / r;
}

// ----------------------------------------------------------------------

double PairPINN::bop_zeta(int it, int jt, int kt, Param *param,
                             double s_ij, double cs_ijk, double r)
{
  double cut;

  cut = ann_fc (r, param);
  return s_ij * small_a[it][jt][kt] * cut
      * square(cs_ijk - small_h[it][jt][kt]);
}

// ----------------------------------------------------------------------

double PairPINN::screen_ijk(int it, int jt, int kt, Param *param, double r)
{
  return  1.0 - ann_fc (r, param) * exp(-lambda[it][jt][kt] * fabs(r));
}

// ----------------------------------------------------------------------

double PairPINN::bop_fr(int itype, int jtype, double r)
{
  // only f_R(r) as in Tersoff

  return exp(big_a[itype][jtype] - alpha[itype][jtype] * r);
}

// ----------------------------------------------------------------------

double PairPINN::bop_fr_d(int itype, int jtype, Param *param, double r)
{
  // derivative of f_R(r)f_c(r) wrt 'r'

  double tmp_exp, tmp_fc, tmp_fc_d;

  tmp_fc = ann_fc (r, param);
  tmp_fc_d = ann_fc_d (r, param);
  tmp_exp = bop_fr (itype, jtype, r);

  return tmp_exp * (tmp_fc_d - tmp_fc * alpha[itype][jtype]) / r;
}


// ----------------------------------------------------------------------

void PairPINN::read_file(char *file)
{
  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = fopen(file, "r");
    if (fp == nullptr) {
      char str[128];
      snprintf(str, 128,"Cannot open PINN potential file %s",file);
      error->one(FLERR,str);
    }
  }

  int nwords, bufsize, nparams;
  char line[MAXLINE],*ptr;
  char **words;
  double *NNparams;

  // first 3 header lines are skipped !!!
  if (comm->me == 0) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  }

  if (comm->me == 0) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    bufsize = strlen(line) + 1;
  }
  MPI_Bcast(&bufsize, 1, MPI_INT, 0, world);
  MPI_Bcast(line, bufsize, MPI_CHAR, 0, world);

  if ((ptr = strchr(line,'#'))) *ptr = '\0';
  nwords = utils::count_words(line);
  if (nwords < 3)
    error->all (FLERR, "Invalid entries");
  sscanf(line, "%d %lf %d", &params->gimethod, &params->giref,
         &params->actfunc);

  if (comm->me == 0) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    bufsize = strlen(line) + 1;
  }
  MPI_Bcast(&bufsize, 1, MPI_INT, 0, world);
  MPI_Bcast(line, bufsize, MPI_CHAR, 0, world);

  sscanf(line, "%d", &nelements);

  elements = memory->create(elements, nelements, 4,"");
  params->mass = memory->create(params->mass, nelements, "");

  // read element name and atomic mass
  for (int i=0; i<nelements; i++) {
    if (comm->me == 0) {
      utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
      bufsize = strlen(line) + 1;
    }
    MPI_Bcast(&bufsize, 1, MPI_INT, 0, world);
    MPI_Bcast(line, bufsize, MPI_CHAR, 0, world);

    nwords = sscanf(line, "%s %lf", elements[i], &params->mass[i]);
    if (nwords != 2) {
      if (comm->me == 0) fclose(fp);
      error->all (FLERR, "Invalid entries");
    }
  }

  if (comm->me == 0) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    bufsize = strlen(line) + 1;
  }
  MPI_Bcast(&bufsize, 1, MPI_INT, 0, world);
  MPI_Bcast(line, bufsize, MPI_CHAR, 0, world);

  if ((ptr = strchr(line,'#'))) *ptr = '\0';
  nwords = utils::count_words (line);
  if (nwords != 5)
    error->all(FLERR, "Too few or many entries in PINN potential file");

  words = new char*[6];
  nwords = 0;
  words[nwords++] = strtok(line," \t\n\r\f");
  while ((words[nwords++] = strtok(nullptr," \t\n\r\f"))) continue;
  params->cut = atof(words[2]);
  params->cut_range = atof(words[3]);
  params->gwidth = atof(words[4]);
  params->cutsq = params->cut * params->cut;

  delete [] words;

  // read number of Legendre-Polynomial orders and the orders.

  if (comm->me == 0) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    bufsize = strlen(line) + 1;
  }
  MPI_Bcast(&bufsize, 1, MPI_INT, 0, world);
  MPI_Bcast(line, bufsize, MPI_CHAR, 0, world);

  if ((ptr = strchr(line,'#'))) *ptr = '\0';
  sscanf(line, "%d", &params->nLPOrders);
  nwords = utils::count_words (line);
  if (nwords != params->nLPOrders + 1) {
    if (comm->me == 0) fclose(fp);
    error->all(FLERR,
               "Incorrect Legendre Polynomial orders in PINN potential file");
  }

  words = new char*[params->nLPOrders + 1];
  nwords = 0;
  strtok(line," \t\n\r\f"); // skip first word
  while ((words[nwords++] = strtok(nullptr," \t\n\r\f"))) continue;

  if (params->nLPOrders > 0) {
    params->LPOrders = memory->create(params->LPOrders, params->nLPOrders, "");
    for (int i=0; i<params->nLPOrders; i++) {
      params->LPOrders[i] = atoi(words[i]);
    }
  } else {
    if (comm->me == 0) fclose(fp);
    error->all(FLERR,
      "Illegal number of Legendre Polynomial orders in PINN potential file");
  }

  delete [] words;

  // read number of positions of gaussian functions and positions

  if (comm->me == 0) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    bufsize = strlen(line) + 1;
  }
  MPI_Bcast(&bufsize, 1, MPI_INT, 0, world);
  MPI_Bcast(line, bufsize, MPI_CHAR, 0, world);

  if ((ptr = strchr(line,'#'))) *ptr = '\0';
  sscanf(line, "%d", &params->nsigmas);
  nwords = utils::count_words (line);
  if (nwords < params->nsigmas + 1) {
    if (comm->me == 0) fclose(fp);
    error->all (FLERR, "Incorrect number of positions in PINN potential file");
  }

  words = new char*[params->nsigmas + 1];
  nwords = 0;
  strtok(line," \t\n\r\f");
  while ((words[nwords++] = strtok(nullptr," \t\n\r\f"))) continue;

  if (params->nsigmas > 0) {
    params->sigmas = memory->create(params->sigmas, params->nsigmas, "");
    for (int i=0; i<params->nsigmas; i++) {
      params->sigmas[i] = atof(words[i]);
    }
  } else {
    if (comm->me == 0) fclose(fp);
    error->all(FLERR, "Illegal number of positions in PINN potential file");
  }

  delete [] words;

  // read number of layers and number of nodes in each layer

  if (comm->me == 0) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    bufsize = strlen(line) + 1;
  }
  MPI_Bcast(&bufsize, 1, MPI_INT, 0, world);
  MPI_Bcast(line, bufsize, MPI_CHAR, 0, world);

  if ((ptr = strchr(line,'#'))) *ptr = '\0';
  sscanf(line, "%d", &params->nlayers);
  nwords = utils::count_words (line);
  if (nwords != params->nlayers + 1)
    error->all (FLERR, "Incorrect layers in PINN potential file");

  words = new char*[params->nlayers + 1];
  nwords = 0;
  strtok(line," \t\n\r\f");
  while ((words[nwords++] = strtok(nullptr," \t\n\r\f"))) continue;

  if (params->nlayers > 0) {
    params->layers = memory->create(params->layers, params->nlayers, "");
    for (int i=0; i<params->nlayers; i++) {
      params->layers[i].nnodes = atoi(words[i]);
    }
  } else {
    if (comm->me == 0) fclose(fp);
    error->all(FLERR,
               "Illegal number of network layers in PINN potential file");
  }

  delete [] words;

  // set flags for PINN Kind 1 and 2

  pinn_kind_flag1 = params->layers[0].nnodes /
      (nelements * (nelements + 1) / 2 * params->nLPOrders * params->nsigmas);

  if (pinn_kind_flag1 == 1) pinn_kind_flag2 = 0;
  else {
    if (pinn_kind_flag1 == nelements) pinn_kind_flag2 = 1;
    else {
      if (comm->me == 0) fclose(fp);
      error->all(FLERR, "Cannot determine PINN Kind");
    }
  }

  // allocate memory for neural network parameters

  for (int i=0; i<params->nlayers; i++) {
    params->layers[i].Weights = nullptr;
    params->layers[i].Biases = nullptr;
    params->layers[i].fdot = nullptr;
    if (i) {
      params->layers[i].Weights =
          memory->create(params->layers[i].Weights,
                         params->layers[i-1].nnodes * params->layers[i].nnodes,
                         "");
      params->layers[i].Biases = memory->create(params->layers[i].Biases,
          params->layers[i].nnodes, "");
      params->layers[i].fdot = memory->create(params->layers[i].fdot,
                                              params->layers[i].nnodes, "");
    }
  }

  // # of baseline BOPs same as # of nodes in the output layer !!!
  nbaseline_bo_params = params->layers[params->nlayers - 1].nnodes;

  // for sanity check
  // first determine the size of the array to store the parameters

  int nbo_params = 4 * nelements;          // As, alphas, Bs and betas
  nbo_params += nelements * nelements + 1; // hs and sigma
  nbo_params += 2 * nelements * nelements; // small a's and lambdas
  nbo_params *= nelements;                 // sets of nelements

  if (nbo_params != nbaseline_bo_params ||
      nbo_params != params->layers[params->nlayers - 1].nnodes) {
    if (comm->me == 0) fclose(fp);
    error->all(FLERR,
           "Illegal number of baseline BO parameters in PINN potential file");
  }

  baseline_bo_params = memory->create(baseline_bo_params,
                                      nbaseline_bo_params, "");

  // read number of baseline BO parameters and the parameters
  // first word in first line must be equal to # of baseline BO parameters
  // followed by parameters spanned over several lines

  if (comm->me == 0) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    bufsize = strlen(line) + 1;
  }

  MPI_Bcast(&bufsize, 1, MPI_INT, 0, world);
  MPI_Bcast(line, bufsize, MPI_CHAR, 0, world);

  if ((ptr = strchr(line,'#'))) *ptr = '\0';
  nwords = utils::count_words (line);

  int i = 0;
  nBase = atoi(strtok(line," \t\n\r\f"));

  while (i < nwords - 1) {
    baseline_bo_params[i++] = atof(strtok(nullptr," \t\n\r\f"));
  }

  // read other lines until all BO parameters are read

  while (i < nbaseline_bo_params) {
    if (comm->me == 0) {
      utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
      bufsize = strlen(line) + 1;
    }

    MPI_Bcast(&bufsize, 1, MPI_INT, 0, world);
    MPI_Bcast(line, bufsize, MPI_CHAR, 0, world);

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = utils::count_words (line);

    for (int j = 0; j < nwords; j++) {
      if (!j) {
        baseline_bo_params[i++] = atof(strtok(line," \t\n\r\f"));
      }
      else {
        baseline_bo_params[i++] = atof(strtok(nullptr," \t\n\r\f"));
      }
    }
  }

  // prepare an array to store weights and biases

  nparams = 0;

  for (i=1; i<params->nlayers; i++) {
    nparams += params->layers[i].nnodes
        + params->layers[i-1].nnodes * params->layers[i].nnodes;
  }

  NNparams = new double[nparams];

  // first read weights and biases in to an array

  i = 0;

  while (i < nparams) {
    if (comm->me == 0) {
      utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
      bufsize = strlen(line) + 1;
    }
    MPI_Bcast(&bufsize, 1, MPI_INT, 0, world);
    MPI_Bcast(line, bufsize, MPI_CHAR, 0, world);

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = utils::count_words (line);

    for (int j = 0; j < nwords; ++j) {
      if (!j) {
        NNparams[i++] = atof(strtok(line," \t\n\r\f"));
      } else {
        NNparams[i++] = atof(strtok(nullptr," \t\n\r\f"));
      }
    }
  }

  if (comm->me == 0) fclose(fp);

  // unpack parameters

  nwords = 0;
  for (int l=1; l<params->nlayers; l++) {
    // weights: m x n matrix
    for (int m=0; m<params->layers[l - 1].nnodes; m++) { // rows
      for (int n=0; n<params->layers[l].nnodes; n++) { // cols
        params->layers[l].Weights[m * params->layers[l].nnodes + n]
            = NNparams[nwords++];
      }
    }
    // biases: vector of size nnodes.
    for (int m=0; m<params->layers[l].nnodes; m++) { // cols
      params->layers[l].Biases[m] = NNparams[nwords++];
    }
  }

  // sanity check !
  if (nwords != nparams)
    error->all(FLERR,"Insufficient weights and biases in PINN potential file");

  if (NNparams) delete [] NNparams;
}

// ----------------------------------------------------------------------

void PairPINN::eval_nnet(double *input, int insize,
                         double *output, int outsize)
{
  // sanity checks !!!
  if (insize != params->layers[0].nnodes ||
      outsize != params->layers[params->nlayers - 1].nnodes)
  {
    error->all(FLERR, "Illegal size(s)");
  }

  int i, j;
  int nrows, ncols;
  double *sum, *prod;

  for (i = 1; i < params->nlayers-1; i++) { // only hidden layers
    if (i == 1) {
      prod = (double *) malloc (params->layers[i].nnodes * sizeof(double));
      sum = (double *) malloc (params->layers[i].nnodes * sizeof(double));
      // Multiply 'input' and 'weights' of the first hidden layer.
      nrows = insize; // # of rows of the matrix which equals size of the vector
      ncols = params->layers[i].nnodes; // # of columns of the matrix
      vec_mult_mat (input, params->layers[i].Weights, prod, nrows, ncols);
    }
    else {
      prod = (double *) realloc(prod, params->layers[i].nnodes
                                * sizeof(double));
      // 'sum' holds result from previous layer !!!
      // Multiply 'sum' and 'weights' of this layer.
      nrows = params->layers[i-1].nnodes; // # of rows of the matrix which equals size of the vector
      ncols = params->layers[i].nnodes; // # of columns of the matrix
      vec_mult_mat (sum, params->layers[i].Weights, prod, nrows, ncols);
      sum = (double *) realloc(sum, params->layers[i].nnodes * sizeof(double));
    }

    // Add 'prod' and 'biases' and save in 'prod'.
    for (int m = 0; m < params->layers[i].nnodes; m++) {
      prod[m] += params->layers[i].Biases[m];
    }

    // Apply activation function to a vector.
    switch (params->actfunc) {
    case 0: // sigmoid function
      for (j=0; j<params->layers[i].nnodes; j++) {
        prod[j] = 1.0 / (1.0 + exp(-prod[j]));
        params->layers[i].fdot[j] = prod[j] * (1.0 - prod[j]);
      }
      break;
    case 1: // 1/2 tanh (x/2)
      for (j=0; j<params->layers[i].nnodes; j++) {
        prod[j] = 1.0 / (1.0 + exp(-prod[j])) - 0.5;
        params->layers[i].fdot[j] = 0.25 - prod[j] * prod[j];
      }
      break;
    case 2: // tanh(x)
      for (j=0; j<params->layers[i].nnodes; j++) {
        prod[j] = tanh(prod[j]);
        params->layers[i].fdot[j] = 1.0 - prod[j] * prod[j];
      }
      break;
    default:
      error->all (FLERR, "Invalid activation function");
    }

    // Copy of 'prod' to 'sum'.
    for (int m = 0; m < params->layers[i].nnodes; m++) {
      sum[m] = prod[m];
    }
  }

  // Output layer does not apply activation function.
  prod = (double *) realloc(prod, params->layers[i].nnodes * sizeof(double));
  nrows = params->layers[i-1].nnodes;
  ncols = params->layers[i].nnodes;
  vec_mult_mat (sum, params->layers[i].Weights, prod, nrows, ncols);

  // Add 'prod' and 'biases' and save in 'output'.
  for (int m = 0; m < params->layers[i].nnodes; m++) {
    output[m] = prod[m] + params->layers[i].Biases[m];
  }

  if (sum) free(sum);
  if (prod) free(prod);
}

void PairPINN::eval_nnet_d(Layer *aNNlayers, double *input,
                             int insize, double *output, int outsize)
{
  // sanity checks !!!
  if (insize != params->layers[0].nnodes ||
      outsize != params->layers[params->nlayers - 1].nnodes)
  {
    error->all(FLERR, "Illegal size(s)");
  }

  int i, j;
  int nrows, ncols;
  double *sum, *prod;

  for (i = 1; i < params->nlayers-1; i++) { // only hidden layers
    if (i == 1) {
      prod = (double *) malloc (params->layers[i].nnodes * sizeof(double));
      sum = (double *) malloc (params->layers[i].nnodes * sizeof(double));
      // Multiply 'input' and 'weights' of this layer.
      nrows = insize;
      ncols = params->layers[i].nnodes;
      vec_mult_mat (input, params->layers[i].Weights, prod, nrows, ncols);
    }
    else {
      prod = (double *) realloc(prod, params->layers[i].nnodes*sizeof(double));
      // 'sum' holds result from previous layer !!!
      nrows = params->layers[i-1].nnodes;
      ncols = params->layers[i].nnodes;
      vec_mult_mat (sum, params->layers[i].Weights, prod, nrows, ncols);
      sum = (double *) realloc(sum, params->layers[i].nnodes * sizeof(double));
    }

    for (j = 0; j < params->layers[i].nnodes; j++) {
      prod[j] = prod[j] * aNNlayers[i].fdot[j];
      sum[j] = prod[j];
    }
  }

  prod = (double *) realloc(prod, params->layers[i].nnodes * sizeof(double));
  nrows = params->layers[i-1].nnodes;
  ncols = params->layers[i].nnodes;
  vec_mult_mat (sum, params->layers[i].Weights, prod, nrows, ncols);
  // Copy 'prod' to 'output'
  for (j = 0; j < outsize; ++j) {
    output[j] = prod[j];
  }

  if (sum) free(sum);
  if (prod) free(prod);
}

void PairPINN::vec_mult_mat(double *vec, double *mat, double *prod,
                            const int nrows, const int ncols)
{
  for (int n = 0; n < ncols; ++n) {
    prod[n] = 0.0;
    for (int m = 0; m < nrows; ++m) {
      prod[n] += vec[m] * mat[n + m * ncols];
    }
  }
}


void PairPINN::compute_atomic_LSP(const int iam, double *gi)
{
  // Computes local structural parameters of an atom.

  double rsqij, rsqik;
  double rij, rik;
  double *Pl;
  double xtmp, ytmp, ztmp;
  double xij, yij, zij;
  double xik, yik, zik;
  int nblock;
  int unit_block;
  int ngi;
  double rc, rc2;
  int j, k, jnum, itype, jtype, ktype;
  int *jlist, *numneigh, **firstneigh;

  // Gis use different cutoff !
  rc = cutmax; // 1.5*cut for PINN !
  rc2 = rc * rc;

  Pl = new double [params->LPOrders[params->nLPOrders - 1] + 1];
  Pl[0] = 1.0;

  // number of blocks of gis of size, (nsigmas * nLPOrders), to fill.
  nblock = nelements * (nelements + 1) / 2;
  unit_block = params->nsigmas * params->nLPOrders;
  ngi = pinn_kind_flag1 * nblock * unit_block;

  double **x = atom->x;
  int *type = atom->type;

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  jnum = numneigh[iam];
  jlist = firstneigh[iam];
  itype = map[type[iam]];

  xtmp = x[iam][0];
  ytmp = x[iam][1];
  ztmp = x[iam][2];

  //for (int n=0; n<params->nsigmas; n++) { // loop over sigmas
  for (int jj = 0; jj < jnum; jj++) {    // loop over sites j
    j = jlist[jj];
    j &= NEIGHMASK;
    jtype = map[type[j]];

    xij = x[j][0] - xtmp;
    yij = x[j][1] - ytmp;
    zij = x[j][2] - ztmp;
    rsqij = xij * xij + yij * yij + zij * zij;

    if (rsqij > rc2) continue;
    rij = sqrt(rsqij);

    for (int kk = 0; kk < jnum; kk++) {  // loop over sites k
      k = jlist[kk];
      k &= NEIGHMASK;
      ktype = map[type[k]];

      xik = x[k][0] - xtmp;
      yik = x[k][1] - ytmp;
      zik = x[k][2] - ztmp;
      rsqik = xik * xik + yik * yik + zik * zik;

      if (rsqik > rc2) continue;
      rik = sqrt(rsqik);

      double costheta = (xij * xik + yij * yik + zij * zik) / (rij * rik);
      Pl[1] = costheta;
      for (int m = 1; m < params->LPOrders[params->nLPOrders - 1]; m++)
        Pl[m + 1] = ((2.0 * m + 1.0) * costheta * Pl[m] - m
                     * Pl[m - 1]) / (m + 1);
      for (int n = 0; n < params->nsigmas; n++) { // loop over sigmas
        double tmp = exp(-square((rij - params->sigmas[n]) / params->gwidth))
            * ann_fc (rij, rc, params->cut_range);
        tmp *= exp(-square((rik - params->sigmas[n]) / params->gwidth))
            * ann_fc (rik, rc, params->cut_range) / square(params->sigmas[n]);
        // loop over Legendre polynomial orders
        for (int p = 0; p < params->nLPOrders; p++) {
          gi[itype * nblock * unit_block * pinn_kind_flag2
              + map_gi[jtype][ktype] * unit_block
              + p * params->nsigmas + n] += Pl[params->LPOrders[p]] * tmp;
        }
      }
    }
  }

  switch(params->gimethod) {
  case 0: // Do nothing.
    break;
  case 2: // Transform all Gis with Inverse Sine Hyperbolic function.
    for (size_t i = 0; i < ngi; i++)
      gi[i] = asinh(gi[i]);
    break;
  default:
    error->all (FLERR, "Invalid gi calculation method");
  }

  if (Pl) delete [] Pl;
}

double PairPINN::ann_fc (double r, Param *par)
{
  double x, y, x4;

  if (r > par->cut) return 0.0;

  x = r - par->cut;
  x = x / par->cut_range;
  x4 = powint(x, 4);
  y = 1.0 + x4;

  return x4 / y;
}

double PairPINN::ann_fc_d (double r, Param *par)
{
  double x,y,x4;
  double x1, y2;

  if (r > par->cut) return 0.0;

  x = r - par->cut;
  x = x / par->cut_range;
  x4 = powint(x, 4);
  y = 1.0 + x4;

  x1 = 4.0 * x * x;
  y2 = y * y;

  return x1 * x / y2 / par->cut_range;
}

double PairPINN::ann_fc (double r, double rc, double hc)
{
  double x, y, x4;

  if (r > rc) return 0.0;

  x = r - rc;
  x = x / hc;
  x4 = powint(x, 4);
  y = 1.0 + x4;

  return x4 / y;
}

// -----------------------------------------------------------------------

double PairPINN::ann_fc_d (double r, double rc, double hc)
{
  double x,y,x4;
  double x1, y2;

  if (r > rc) return 0.0;

  x = r - rc;
  x = x / hc;
  x4 = powint(x, 4);
  y = 1.0 + x4;

  x1 = 4.0 * x * x;
  y2 = y * y;

  return x1 * x / y2 / hc;
}

// -----------------------------------------------------------------------

double PairPINN::ann_fs(double r, int s, Param *par)
{
  double tmpexp = exp(-square((r - par->sigmas[s]) / par->gwidth))
      / par->sigmas[s];
  double tmp_fc = ann_fc(r, par);

  return tmpexp * tmp_fc;
}

// -----------------------------------------------------------------------

double PairPINN::pinn_fs(const double r, const int s)
{
  double tmpexp = exp(-square((r - params->sigmas[s]) / params->gwidth))
      / params->sigmas[s];
  double tmp_fc = ann_fc(r, cutmax, params->cut_range);

  return tmpexp * tmp_fc;
}

// -----------------------------------------------------------------------

double PairPINN::ann_fs_d(double r, int s, Param *par)
{
  double tmp_exp = exp(-square((r - par->sigmas[s]) / par->gwidth))
      / par->sigmas[s];
  double tmp1 = tmp_exp * ann_fc_d(r, par);
  double tmp2 = -2.0 * tmp_exp * ann_fc(r, par) * (r - par->sigmas[s])
      / par->gwidth / par->gwidth;

  return (tmp1 + tmp2);
}

// -----------------------------------------------------------------------

double PairPINN::pinn_fs_d(const double r, const int s)
{
  double tmp_exp = exp(-square((r - params->sigmas[s]) / params->gwidth))
      / params->sigmas[s];
  double tmp1 = tmp_exp * ann_fc_d(r, cutmax, params->cut_range);
  double tmp2 = -2.0 * tmp_exp * ann_fc(r, cutmax, params->cut_range) *
      (r - params->sigmas[s]) / params->gwidth / params->gwidth;

  return (tmp1 + tmp2);
}

// -----------------------------------------------------------------------

void PairPINN::unpack_vec_to_bop_sets(const double *pvec, const int nsize)
{
  // ------------------------------------------
  // Given parameter vector is unpacked
  // into corresponding sets of straight
  // BO parameters in a specified order.
  // The order is important !!!
  //
  // Current order:
  //
  // A, alpha, B, beta, h, sigma, a, lambda
  // ------------------------------------------


  // check size of input vector
  if (nsize <= 0) error->all (FLERR, "Invalid size (n)");

  // check if multidimensional arrays are allocated
  if (!big_a && !lambda) error->all (FLERR, "Array(s) not allocated");

  int pcount = 0;

  for (int n = 0; n < nelements; n++) { // loop over number of types
    // As
    for (int i = 0; i < nelements; i++) {
      big_a[n][i] = pvec[pcount++];
    }

    // alphas
    for (int i = 0; i < nelements; i++) {
      alpha[n][i] = pvec[pcount++];
    }

    // Bs
    for (int i = 0; i < nelements; i++) {
      big_b[n][i] = pvec[pcount++];
    }

    // betas
    for (int i = 0; i < nelements; i++) {
      beta[n][i] = pvec[pcount++];
    }

    // hs
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        small_h[n][i][j] = pvec[pcount++];
      }
    }

    // sigmas
    sigma[n] = pvec[pcount++];

    // a's
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        small_a[n][i][j] = pvec[pcount++];
      }
    }

    // lambdas
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        lambda[n][i][j] = pvec[pcount++];
      }
    }
  }

  if (pcount != nsize) {
    error->all (FLERR, "Insufficient number of BO parameters");
  }
}

void PairPINN::iderivs_wrt_bops(const int myID, const double ep)
{
  double xtmp, ytmp, ztmp;
  int j, k, jnum, itype, jtype, ktype;
  int *jlist, *numneigh, **firstneigh;
  double r, rsqij, rsqik, rij, rik, rjk;
  double drij[3], drik[3], drjk[3], dril[3], drkl[3];

  double **x = atom->x;
  int *type = atom->type;

  // initialize arrays to store derivatives w.r.t. BO parameters

  memset(deriv_wrt_big_a[0], 0, nelements * nelements * sizeof (double));
  memset(deriv_wrt_alpha[0], 0, nelements * nelements * sizeof (double));
  memset(deriv_wrt_big_b[0], 0, nelements * nelements * sizeof (double));
  memset(deriv_wrt_beta[0], 0, nelements * nelements * sizeof (double));
  memset(deriv_wrt_small_h[0][0], 0, nelements * nelements * nelements
      * sizeof (double));
  memset(deriv_wrt_sigma, 0, nelements * sizeof (double));
  memset(deriv_wrt_small_a[0][0], 0, nelements * nelements * nelements
      * sizeof (double));
  memset(deriv_wrt_lambda[0][0], 0, nelements * nelements * nelements
      * sizeof (double));

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  jnum = numneigh[myID];
  jlist = firstneigh[myID];
  itype = map[type[myID]];

  xtmp = x[myID][0];
  ytmp = x[myID][1];
  ztmp = x[myID][2];

  // compute derivatives of E_i w.r.t. local BO parameters

  for (int jj = 0; jj < jnum; jj++) { // loop over neighbors of current atom
    j = jlist[jj];
    j &= NEIGHMASK;
    jtype = map[type[j]];

    drij[0] = x[j][0] - xtmp;
    drij[1] = x[j][1] - ytmp;
    drij[2] = x[j][2] - ztmp;

    rsqij = dot_product (drij, drij, 3);

    if (rsqij > params->cutsq) continue;
    rij = sqrt(rsqij);

    double tmp = 0.5 * bop_fr (itype, jtype, rij) * ann_fc(rij, params);
    deriv_wrt_big_a[itype][jtype] += tmp;
    deriv_wrt_alpha[itype][jtype] += -rij * tmp;

    tmp = 0.5 * Sij[jj] * bij[jj] * bop_fa (itype, jtype, rij)
        * ann_fc (rij, params);
    deriv_wrt_big_b[itype][jtype] += tmp;
    deriv_wrt_beta[itype][jtype] += -rij * tmp;

    double dbij = powint(bij[jj], 3) * copysign(1.0, 1.0 + Zij[jj]);

    double prefactor_dzij_1 = -0.25 * Sij[jj] * dbij *
        bop_fa (itype, jtype, rij) * ann_fc (rij, params);

    double prefactor_dzij_2 = 0.25 * sigma[itype] / ep * Sij[jj] *
        dbij * ann_fc (rij, params);

    double prefactor_dsij_dli_1 = 0.5 * bij[jj] * bop_fa (itype, jtype, rij)
        * ann_fc (rij, params);
    double prefactor_dsij_dli_2 = -0.5 * sigma[itype] / ep * bij[jj]
        * ann_fc (rij, params);

    for (int kk = 0; kk < jnum; ++kk) {
      if (jj == kk) continue;
      k = jlist[kk];
      k &= NEIGHMASK;
      ktype = map[type[k]];

      drik[0] = x[k][0] - xtmp;
      drik[1] = x[k][1] - ytmp;
      drik[2] = x[k][2] - ztmp;

      rsqik = dot_product (drik, drik, 3);
      rik = sqrt(rsqik);

      drjk[0] = drik[0] - drij[0]; // rj ----> rk
      drjk[1] = drik[1] - drij[1];
      drjk[2] = drik[2] - drij[2];

      rjk = sqrt(dot_product (drjk, drjk, 3));

      r = rik + rjk - rij;

      if (r < params->cut) {
        deriv_wrt_lambda[itype][jtype][ktype] +=
            (prefactor_dsij_dli_1 + prefactor_dsij_dli_2)
            * Sij[jj] / Sijk[jj][kk] * r *exp(-lambda[itype][jtype][ktype] * r)
            * ann_fc (r, params);
      }

      if (rsqik < params->cutsq) {
        deriv_wrt_small_a[itype][jtype][ktype] +=
            (prefactor_dzij_1 + prefactor_dzij_2) * Sij[kk]
            * square(CSijk[jj][kk] - small_h[itype][jtype][ktype])
            * ann_fc(rik, params);
        deriv_wrt_small_h[itype][jtype][ktype] +=
            -(prefactor_dzij_1 + prefactor_dzij_2) *
            2.0 * Sij[kk] * small_a[itype][jtype][ktype] *
            (CSijk[jj][kk] - small_h[itype][jtype][ktype])
            * ann_fc(rik, params);

        double prefactor_dsik_dli_1 = -0.25 * Sij[jj]
            * bop_fa(itype, jtype, rij) * ann_fc (rij, params) * dbij
            * small_a[itype][jtype][ktype]
            * square(CSijk[jj][kk] - small_h[itype][jtype][ktype])
            * ann_fc (rik, params);

        double prefactor_dsik_dli_2 = 0.25 * sigma[itype] / ep * Sij[jj] * dbij
            * ann_fc (rij, params) * small_a[itype][jtype][ktype]
            * square(CSijk[jj][kk] - small_h[itype][jtype][ktype])
            * ann_fc (rik, params);

        for (int ll = 0; ll < jnum; ++ll) {
          if (ll == kk) continue;
          int l = jlist[ll];
          l &= NEIGHMASK;
          int ltype = map[type[l]];

          dril[0] = x[l][0] - xtmp;
          dril[1] = x[l][1] - ytmp;
          dril[2] = x[l][2] - ztmp;

          double ril = sqrt(dot_product (dril, dril, 3));

          drkl[0] = dril[0] - drik[0];
          drkl[1] = dril[1] - drik[1];
          drkl[2] = dril[2] - drik[2];

          double rkl = sqrt(dot_product (drkl, drkl, 3));

          r = ril + rkl - rik;

          if (r < params->cut) {
            deriv_wrt_lambda[itype][ktype][ltype] +=
                (prefactor_dsik_dli_1 + prefactor_dsik_dli_2)
                * Sij[kk] / Sijk[kk][ll] * r
                * exp(-lambda[itype][ktype][ltype] * r)
                * ann_fc (r, params);
          }
        }
      }
    }
  }

  deriv_wrt_sigma[itype] += -ep;
}

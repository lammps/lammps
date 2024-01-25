// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Environment Dependent Interatomic Potential

   Contributing author: Chao Jiang
------------------------------------------------------------------------- */

#include "pair_edip_multi.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace MathExtra;

static constexpr int MAXLINE = 1024;
static constexpr int DELTA = 4;

static const char cite_pair_edip[] =
  "pair edip/multi: doi:10.1103/PhysRevB.86.144118, doi:10.1088/0953-8984/22/3/035802\n\n"
  "@article{cjiang2012\n"
  " author    = {Jian, Chao and Morgan, Dane, and Szlufarska, Izabella},\n"
  " title     = {Carbon Tri-Interstitial Defect: {A} Model for {D$_{\\mathrm{II}}$} Center},\n"
  " journal   = {Phys.\\ Rev.~B},\n"
  " volume    = {86},\n"
  " pages     = {144118},\n"
  " year      = {2012},\n"
  "}\n\n"
  "@article{lpizzagalli2010,\n"
  " author    = {G. Lucas and M. Bertolus and L. Pizzagalli},\n"
  " journal   = {J.~Phys.\\ Condens.\\ Matter},\n"
  " volume    = {22},\n"
  " number    = 3,\n"
  " pages     = {035802},\n"
  " year      = {2010},\n"
  "}\n\n";

// max number of interaction per atom for f(Z) environment potential

static constexpr int leadDimInteractionList = 64;

static inline void costheta_d(const double *dr_ij, const double r_ij,
                              const double *dr_ik, const double r_ik,
                              double *dri, double *drj, double *drk)
{
  const double costheta = dot3(dr_ij, dr_ik) / r_ij / r_ik;
  scaleadd3(1 / r_ij / r_ik, dr_ik, -costheta / r_ij / r_ij, dr_ij, drj);
  scaleadd3(1 / r_ij / r_ik, dr_ij, -costheta / r_ik / r_ik, dr_ik, drk);
  scaleadd3(-1, drj, -1, drk, dri);
}

/* ---------------------------------------------------------------------- */

PairEDIPMulti::PairEDIPMulti(LAMMPS *lmp) : Pair(lmp), preForceCoord(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_edip);

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;

  params = nullptr;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairEDIPMulti::~PairEDIPMulti()
{
  memory->destroy(params);
  memory->destroy(elem3param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
  deallocatePreLoops();
}

/* ---------------------------------------------------------------------- */

void PairEDIPMulti::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  double xtmp,ytmp,ztmp,evdwl;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int preForceCoord_counter;

  double zeta_i;
  double dzetair;
  double fpair;
  double costheta;
  double dpairZ,dtripleZ;

  // eflag != 0 means compute energy contributions in this step
  // vflag != 0 means compute virial contributions in this step

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;//total number of atoms in the cell
  ilist = list->ilist;//list of atoms
  numneigh = list->numneigh;//number of near neighbors
  firstneigh = list->firstneigh;//list of neighbors

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    zeta_i = 0.0;
    int numForceCoordPairs = 0;

    i = ilist[ii];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // all the neighbors of atom i

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // pre-loop to compute environment coordination f(Z)

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      double delx, dely, delz, r_ij;

      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      r_ij = delx * delx + dely * dely + delz * delz;

      jtype = map[type[j]];
      const Param &param = params[elem3param[itype][jtype][jtype]];

      if (r_ij > param.cutsq) continue;

      r_ij = sqrt(r_ij);

      // zeta and its derivative dZ/dr

      if (r_ij < param.cutoffC) zeta_i += 1.0;
      else {
        double f, fdr;
        edip_fc(r_ij, param, f, fdr);
        zeta_i += f;
        dzetair = -fdr / r_ij;

        preForceCoord_counter=numForceCoordPairs*5;
        preForceCoord[preForceCoord_counter+0]=dzetair;
        preForceCoord[preForceCoord_counter+1]=delx;
        preForceCoord[preForceCoord_counter+2]=dely;
        preForceCoord[preForceCoord_counter+3]=delz;
        preForceCoord[preForceCoord_counter+4]=j;
        numForceCoordPairs++;
      }
    }

    // two-body interactions

    dpairZ=0;
    dtripleZ=0;

    for (jj = 0; jj < jnum; jj++) {
      double dr_ij[3], r_ij, f_ij[3];

      j = jlist[jj];
      j &= NEIGHMASK;

      dr_ij[0] = x[j][0] - xtmp;
      dr_ij[1] = x[j][1] - ytmp;
      dr_ij[2] = x[j][2] - ztmp;
      r_ij = dr_ij[0]*dr_ij[0] + dr_ij[1]*dr_ij[1] + dr_ij[2]*dr_ij[2];

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (r_ij > params[ijparam].cutsq) continue;

      r_ij = sqrt(r_ij);

      // potential energy and force
      // since pair i-j is different from pair j-i, double counting is
      // already considered in constructing the potential

      double fdr, fdZ;
      edip_pair(r_ij, zeta_i, params[ijparam], evdwl, fdr, fdZ);
      fpair = -fdr / r_ij;
      dpairZ += fdZ;

      f[i][0] -= fpair * dr_ij[0];
      f[i][1] -= fpair * dr_ij[1];
      f[i][2] -= fpair * dr_ij[2];

      f[j][0] += fpair * dr_ij[0];
      f[j][1] += fpair * dr_ij[1];
      f[j][2] += fpair * dr_ij[2];

      if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, -dr_ij[0], -dr_ij[1], -dr_ij[2]);

      // three-body Forces

      for (kk = jj + 1; kk < jnum; kk++) {
        double dr_ik[3], r_ik, f_ik[3];

        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = x[k][0] - xtmp;
        dr_ik[1] = x[k][1] - ytmp;
        dr_ik[2] = x[k][2] - ztmp;
        r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        r_ik = sqrt(r_ik);

        costheta=dot3(dr_ij, dr_ik) / r_ij / r_ik;

        double v1, v2, v3, v4, v5, v6, v7;

        edip_fcut3(r_ij, params[ijparam], v1, v2);
        edip_fcut3(r_ik, params[ikparam], v3, v4);
        edip_h(costheta, zeta_i, params[ijkparam], v5, v6, v7);

        // potential energy and forces
        evdwl = v1 * v3 * v5;
        dtripleZ += v1 * v3 * v7;

        double dri[3], drj[3], drk[3];
        double dhl, dfr;

        dhl = v1 * v3 * v6;

        costheta_d(dr_ij, r_ij, dr_ik, r_ik, dri, drj, drk);

        f_ij[0] = -dhl * drj[0];
        f_ij[1] = -dhl * drj[1];
        f_ij[2] = -dhl * drj[2];
        f_ik[0] = -dhl * drk[0];
        f_ik[1] = -dhl * drk[1];
        f_ik[2] = -dhl * drk[2];

        dfr = v2 * v3 * v5;
        fpair = -dfr / r_ij;

        f_ij[0] += fpair * dr_ij[0];
        f_ij[1] += fpair * dr_ij[1];
        f_ij[2] += fpair * dr_ij[2];

        dfr = v1 * v4 * v5;
        fpair = -dfr / r_ik;

        f_ik[0] += fpair * dr_ik[0];
        f_ik[1] += fpair * dr_ik[1];
        f_ik[2] += fpair * dr_ik[2];

        f[j][0] += f_ij[0];
        f[j][1] += f_ij[1];
        f[j][2] += f_ij[2];

        f[k][0] += f_ik[0];
        f[k][1] += f_ik[1];
        f[k][2] += f_ik[2];

        f[i][0] -= f_ij[0] + f_ik[0];
        f[i][1] -= f_ij[1] + f_ik[1];
        f[i][2] -= f_ij[2] + f_ik[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,f_ij,f_ik,dr_ij,dr_ik);
      }
    }

    // forces due to environment coordination f(Z)
    for (int idx = 0; idx < numForceCoordPairs; idx++) {
      double delx, dely, delz;

      preForceCoord_counter = idx * 5;
      dzetair = preForceCoord[preForceCoord_counter+0];
      delx = preForceCoord[preForceCoord_counter+1];
      dely = preForceCoord[preForceCoord_counter+2];
      delz = preForceCoord[preForceCoord_counter+3];
      j = static_cast<int> (preForceCoord[preForceCoord_counter+4]);

      dzetair *= (dpairZ + dtripleZ);

      f[j][0] += dzetair * delx;
      f[j][1] += dzetair * dely;
      f[j][2] += dzetair * delz;

      f[i][0] -= dzetair * delx;
      f[i][1] -= dzetair * dely;
      f[i][2] -= dzetair * delz;

      evdwl = 0.0;
      if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, dzetair, -delx, -dely, -delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

double sqr(double x)
{
  return x * x;
}

//pair Vij, partial derivatives dVij(r,Z)/dr and dVij(r,Z)/dZ
void PairEDIPMulti::edip_pair(double r, double z, const Param &param, double &eng,
                         double &fdr, double &fZ)
{
  double A = param.A;
  double B = param.B;
  double rho = param.rho;
  double beta = param.beta;
  double v1,v2,v3,v4;

  v1 = pow(B / r, rho);
  v2 = exp(-beta * z * z);
  edip_fcut2(r, param, v3, v4);

  eng = A * (v1 - v2) * v3;
  fdr = A * (v1 - v2) * v4 + A * (-rho * v1 / r) * v3;
  fZ = A * (2 * beta * z * v2) * v3;
}

//function fc(r) in calculating coordination Z and derivative fc'(r)
void PairEDIPMulti::edip_fc(double r, const Param &param, double &f, double &fdr)
{
  double a = param.cutoffA;
  double c = param.cutoffC;
  double alpha = param.alpha;
  double x;
  double v1, v2;

  if (r < c + 1E-6)
  {
    f=1.0;
    fdr=0.0;
    return;
  }

  if (r > a - 1E-6)
  {
    f=0.0;
    fdr=0.0;
    return;
  }

  x = (a - c) / (r - c);
  v1 = x * x * x;
  v2 = 1.0 / (1.0 - v1);

  f = exp(alpha * v2);
  fdr = (3.0 * x * v1 / (a - c)) * (-alpha * v2 * v2) * f;
}

//cut-off function for Vij and its derivative fcut2'(r)
void PairEDIPMulti::edip_fcut2(double r, const Param &param, double &f, double &fdr)
{
  double sigma = param.sigma;
  double a = param.cutoffA;
  double v1;

  if (r > a - 1E-6)
  {
    f=0.0;
    fdr=0.0;
    return;
  }

  v1 = 1.0 / (r - a);
  f = exp(sigma * v1);
  fdr = (-sigma * v1 * v1) * f;
}

//function tau(Z) and its derivative tau'(Z)
void PairEDIPMulti::edip_tau(double z, const Param &param, double &f, double &fdZ)
{
  double u1 = param.u1;
  double u2 = param.u2;
  double u3 = param.u3;
  double u4 = param.u4;
  double v1, v2;

  v1 = exp(-u4 * z);
  v2 = exp(-2.0 * u4 * z);

  f = u1 + u2 * u3 * v1 - u2 * v2;
  fdZ = -u2 * u3 * u4 * v1 + 2.0 * u2 * u4 * v2;
}

//function h(l,Z) and its partial derivatives dh(l,Z)/dl and dh(l,Z)/dZ
void PairEDIPMulti::edip_h(double l, double z, const Param &param, double &f,
                      double &fdl, double &fdZ)
{
  double lambda = param.lambda;
  double eta = param.eta;
  double Q0 = param.Q0;
  double mu = param.mu;
  double Q, QdZ, Tau, TaudZ;
  double u2, du2l, du2Z;
  double v1, v2, v3;

  //function Q(Z)
  Q = Q0 * exp(-mu * z);
  //derivative Q'(Z)
  QdZ= -mu * Q;

  edip_tau(z, param, Tau, TaudZ);

  v1 = sqr(l + Tau);
  u2 = Q * v1;
  v2 = exp(-u2);

  f = lambda * (1 - v2 + eta * u2);

  //df/du2
  v3 = lambda * (v2 + eta);

  //du2/dl
  du2l = Q * 2 * (l + Tau);
  fdl = v3 * du2l;

  //du2/dZ
  du2Z = QdZ * v1 + Q * 2 * (l + Tau) * TaudZ;
  fdZ = v3 * du2Z;
}

//cut-off function for Vijk and its derivative fcut3'(r)
void PairEDIPMulti::edip_fcut3(double r, const Param &param, double &f, double &fdr)
{
  double gamma = param.gamma;
  double a = param.cutoffA;
  double v1;

  if (r > a - 1E-6)
  {
    f=0.0;
    fdr=0.0;
    return;
  }

  v1 = 1.0 / (r - a);
  f = exp(gamma * v1);
  fdr = (-gamma * v1 * v1) * f;
}

/* ----------------------------------------------------------------------
   pre-calculated structures
------------------------------------------------------------------------- */

void PairEDIPMulti::allocatePreLoops()
{
  int nthreads = comm->nthreads;

  memory->create(preForceCoord,5*nthreads*leadDimInteractionList,"edip:preForceCoord");
}

/* ----------------------------------------------------------------------
   deallocate preLoops
------------------------------------------------------------------------- */

void PairEDIPMulti::deallocatePreLoops()
{
  memory->destroy(preForceCoord);
}

/* ---------------------------------------------------------------------- */

void PairEDIPMulti::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEDIPMulti::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairEDIPMulti::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg-3,arg+3);

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup();

  // (re-)allocate tables and internal structures

  deallocatePreLoops();
  allocatePreLoops();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEDIPMulti::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style edip/multi requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style edip/multi requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEDIPMulti::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairEDIPMulti::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "edip");
    char *line;

    while ((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer values(line);
        std::string iname = values.next_string();
        std::string jname = values.next_string();
        std::string kname = values.next_string();

        // ielement,jelement,kelement = 1st args
        // if all 3 args are in element list, then parse this line
        // else skip to next entry in file
        int ielement, jelement, kelement;

        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;
        if (ielement == nelements) continue;
        for (jelement = 0; jelement < nelements; jelement++)
          if (jname == elements[jelement]) break;
        if (jelement == nelements) continue;
        for (kelement = 0; kelement < nelements; kelement++)
          if (kname == elements[kelement]) break;
        if (kelement == nelements) continue;

        // load up parameter settings and error check their values

        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                              "pair:params");

          // make certain all addional allocated storage is initialized
          // to avoid false positives when checking with valgrind

          memset(params + nparams, 0, DELTA*sizeof(Param));
        }

        params[nparams].ielement = ielement;
        params[nparams].jelement = jelement;
        params[nparams].kelement = kelement;
        params[nparams].A = values.next_double();
        params[nparams].B = values.next_double();
        params[nparams].cutoffA = values.next_double();
        params[nparams].cutoffC = values.next_double();
        params[nparams].alpha = values.next_double();
        params[nparams].beta = values.next_double();
        params[nparams].eta = values.next_double();
        params[nparams].gamma = values.next_double();
        params[nparams].lambda = values.next_double();
        params[nparams].mu = values.next_double();
        params[nparams].rho = values.next_double();
        params[nparams].sigma = values.next_double();
        params[nparams].Q0 = values.next_double();
        params[nparams].u1 = values.next_double();
        params[nparams].u2 = values.next_double();
        params[nparams].u3 = values.next_double();
        params[nparams].u4 = values.next_double();
      } catch (std::exception &e) {
        error->one(FLERR, "Error reading EDIP potential file: {}", e.what());
      }

      if (params[nparams].A < 0.0 || params[nparams].B < 0.0 || params[nparams].cutoffA < 0.0 ||
          params[nparams].cutoffC < 0.0 || params[nparams].alpha < 0.0 ||
          params[nparams].beta < 0.0 || params[nparams].eta < 0.0 || params[nparams].gamma < 0.0 ||
          params[nparams].lambda < 0.0 || params[nparams].mu < 0.0 || params[nparams].rho < 0.0 ||
          params[nparams].sigma < 0.0)
        error->all(FLERR, "Illegal EDIP parameter");

      nparams++;
    }
  }
  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if (comm->me != 0)
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");

  MPI_Bcast(params, maxparam*sizeof(Param), MPI_BYTE, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairEDIPMulti::setup()
{
  int i,j,k,m,n;
  double rtmp;

  // set elem3param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem3param);
  memory->create(elem3param,nelements,nelements,nelements,"pair:elem3param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has a duplicate entry for: {} {} {}",
                                   elements[i], elements[j], elements[k]);
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry for: {} {} {}",
                              elements[i], elements[j], elements[k]);
        elem3param[i][j][k] = n;
      }

  // set cutoff square

  for (m = 0; m < nparams; m++) {
    params[m].cutsq = params[m].cutoffA*params[m].cutoffA;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cutsq);
    if (rtmp > cutmax) cutmax = rtmp;
  }
}

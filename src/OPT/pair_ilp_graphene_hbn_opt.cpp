/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Optimization author1: Ping Gao (National Supercomputing Center in Wuxi, China)
   e-mail: qdgaoping at gmail dot com
   Optimization author2: Xiaohui Duan (National Supercomputing Center in Wuxi, China)
   e-mail: sunrise_duan at 126 dot com

   Provides some bugfixes and performance optimizations in this potential.
*/

#include "pair_ilp_graphene_hbn_opt.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "interlayer_taper.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace InterLayer;

static const char cite_ilp_cur[] =
    "ilp/graphene/hbn/opt potential doi:10.1145/3458817.3476137\n"
    "@inproceedings{gao2021lmff\n"
    " author = {Gao, Ping and Duan, Xiaohui and Others},\n"
    " title = {LMFF: Efficient and Scalable Layered Materials Force Field on Heterogeneous "
    "Many-Core Processors},\n"
    " year = {2021},\n"
    " isbn = {9781450384421},\n"
    " publisher = {Association for Computing Machinery},\n"
    " address = {New York, NY, USA},\n"
    " url = {https://doi.org/10.1145/3458817.3476137},\n"
    " doi = {10.1145/3458817.3476137},\n"
    " booktitle = {Proceedings of the International Conference for High Performance Computing, "
    "Networking, Storage and Analysis},\n"
    " articleno = {42},\n"
    " numpages = {14},\n"
    " location = {St. Louis, Missouri},\n"
    " series = {SC'21},\n"
    "}\n\n";

static bool check_vdw(tagint itag, tagint jtag, double *xi, double *xj);

/* ---------------------------------------------------------------------- */

PairILPGrapheneHBNOpt::PairILPGrapheneHBNOpt(LAMMPS *lmp) :
    PairILPGrapheneHBN(lmp), layered_neigh(nullptr), first_layered_neigh(nullptr),
    num_intra(nullptr), num_inter(nullptr), num_vdw(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_ilp_cur);

  single_enable = 0;
  inum_max = 0;
  jnum_max = 0;
}

/* ---------------------------------------------------------------------- */

PairILPGrapheneHBNOpt::~PairILPGrapheneHBNOpt()
{
  memory->destroy(layered_neigh);
  memory->sfree(first_layered_neigh);
  memory->destroy(num_intra);
  memory->destroy(num_inter);
  memory->destroy(num_vdw);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairILPGrapheneHBNOpt::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style ilp/graphene/hbn requires newton pair on");
  if (!atom->molecule_flag)
    error->all(FLERR, "Pair style ilp/graphene/hbn requires atom attribute molecule");

  // It seems that ghost neighbors is required for some historical reason, and it is not needed now

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

void PairILPGrapheneHBNOpt::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  pvector[0] = pvector[1] = 0.0;

  if (neighbor->ago == 0) update_internal_list();

  if (eflag_global || eflag_atom) {
    if (vflag_either) {
      if (tap_flag) {
        eval<1, 1, 1>();
      } else {
        eval<1, 1, 0>();
      }
    } else {
      if (tap_flag) {
        eval<1, 0, 1>();
      } else {
        eval<1, 0, 0>();
      }
    }
  } else {
    if (vflag_either) {
      if (tap_flag) {
        eval<0, 1, 1>();
      } else {
        eval<0, 1, 0>();
      }
    } else {
      if (tap_flag) {
        eval<0, 0, 1>();
      } else {
        eval<0, 0, 0>();
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

template <int EFLAG, int VFLAG_EITHER, int TAP_FLAG> void PairILPGrapheneHBNOpt::eval()
{
  constexpr int EVFLAG = EFLAG || VFLAG_EITHER;
  int i, j, ii, jj, inum, itype, itype_map, jtype, k, kk;
  double prodnorm1, fkcx, fkcy, fkcz;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair, fpair1;
  double rsq, r, Rcut, rhosq1, exp0, exp1, Tap, dTap, Vilp;
  double frho1, Erep, fsum, rdsq1;
  int *ilist;
  evdwl = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double fp1[3] = {0.0, 0.0, 0.0};
  double delki[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};
  double dproddni[3] = {0.0, 0.0, 0.0};
  double cij;

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    itype_map = map[type[i]];
    int *jlist_intra = first_layered_neigh[i];
    int *jlist_inter = first_layered_neigh[i] + num_intra[i];
    int jnum_intra = num_intra[i];
    int jnum_inter = num_inter[i];
    int jnum_vdw = num_vdw[i];
    int ILP_neigh[3];
    int ILP_nneigh = 0;
    for (jj = 0; jj < jnum_intra; jj++) {
      j = jlist_intra[jj];

      jtype = map[type[j]];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutILPsq[itype_map][jtype]) {
        if (ILP_nneigh >= 3)
          error->one(FLERR, "There are too many neighbors for calculating normals");
        ILP_neigh[ILP_nneigh++] = j;
      }
    }    // loop over jj

    dproddni[0] = 0.0;
    dproddni[1] = 0.0;
    dproddni[2] = 0.0;

    double norm[3], dnormdxi[3][3], dnormdxk[3][3][3];
    calc_single_normal(i, ILP_neigh, ILP_nneigh, norm, dnormdxi, dnormdxk);

    for (jj = 0; jj < jnum_inter; jj++) {
      j = jlist_inter[jj];
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      // only include the interaction between different layers
      if (rsq < cutsq[itype][jtype]) {

        int iparam_ij = elem2param[map[itype]][map[jtype]];
        Param &p = params[iparam_ij];

        r = sqrt(rsq);
        double r2inv = 1.0 / rsq;
        double rinv = r * r2inv;
        // turn on/off taper function
        if (TAP_FLAG) {
          Rcut = sqrt(cutsq[itype][jtype]);
          Tap = calc_Tap(r, Rcut);
          dTap = calc_dTap(r, Rcut);
        } else {
          Tap = 1.0;
          dTap = 0.0;
        }

        // Calculate the transverse distance
        prodnorm1 = norm[0] * delx + norm[1] * dely + norm[2] * delz;
        rhosq1 = rsq - prodnorm1 * prodnorm1;    // rho_ij
        rdsq1 = rhosq1 * p.delta2inv;            // (rho_ij/delta)^2

        // store exponents
        exp0 = exp(-p.lambda * (r - p.z0));
        exp1 = exp(-rdsq1);

        frho1 = exp1 * p.C;
        Erep = 0.5 * p.epsilon + frho1;

        Vilp = exp0 * Erep;

        // derivatives
        fpair = p.lambda * exp0 * rinv * Erep;
        fpair1 = 2.0 * exp0 * frho1 * p.delta2inv;
        fsum = fpair + fpair1;
        // derivatives of the product of rij and ni, the result is a vector

        fp1[0] = prodnorm1 * norm[0] * fpair1;
        fp1[1] = prodnorm1 * norm[1] * fpair1;
        fp1[2] = prodnorm1 * norm[2] * fpair1;

        fkcx = (delx * fsum - fp1[0]) * Tap - Vilp * dTap * delx * rinv;
        fkcy = (dely * fsum - fp1[1]) * Tap - Vilp * dTap * dely * rinv;
        fkcz = (delz * fsum - fp1[2]) * Tap - Vilp * dTap * delz * rinv;

        f[i][0] += fkcx;
        f[i][1] += fkcy;
        f[i][2] += fkcz;
        f[j][0] -= fkcx;
        f[j][1] -= fkcy;
        f[j][2] -= fkcz;

        cij = -prodnorm1 * fpair1 * Tap;
        dproddni[0] += cij * delx;
        dproddni[1] += cij * dely;
        dproddni[2] += cij * delz;

        if (EFLAG) pvector[1] += evdwl = Tap * Vilp;
        if (EVFLAG)
          ev_tally_xyz(i, j, nlocal, newton_pair, evdwl, 0.0, fkcx, fkcy, fkcz, delx, dely, delz);

        /* ----------------------------------------------------------------------
	   van der Waals forces and energy
	   ------------------------------------------------------------------------- */
        if (jj >= jnum_vdw) continue;
        double r6inv = r2inv * r2inv * r2inv;
        double r8inv = r6inv * r2inv;

        //double TSvdw = 1.0 + exp(-p.d * (r * p.seffinv - 1.0));
        double TSvdw = 1.0 + exp(-p.d * (r * 1.0 / p.seff - 1.0));
        double TSvdwinv = 1.0 / TSvdw;
        double TSvdw2inv = TSvdwinv * TSvdwinv;    //pow(TSvdw, -2.0);
        Vilp = -p.C6 * r6inv * TSvdwinv;

        fpair = -6.0 * p.C6 * r8inv * TSvdwinv +
            p.C6 * p.d * 1.0 / p.seff * (TSvdw - 1.0) * TSvdw2inv * r8inv * r;
        //p.C6 * p.d * p.seffinv * (TSvdw - 1.0) * TSvdw2inv * r8inv * r;
        fsum = fpair * Tap - Vilp * dTap * rinv;

        double fvdwx = fsum * delx;
        double fvdwy = fsum * dely;
        double fvdwz = fsum * delz;

        f[i][0] += fvdwx;
        f[i][1] += fvdwy;
        f[i][2] += fvdwz;
        f[j][0] -= fvdwx;
        f[j][1] -= fvdwy;
        f[j][2] -= fvdwz;

        if (EFLAG) pvector[0] += evdwl = Tap * Vilp;
        if (EVFLAG)
          ev_tally_xyz(i, j, nlocal, newton_pair, evdwl, 0.0, fvdwx, fvdwy, fvdwz, delx, dely,
                       delz);
      }
    }    // loop over jj

    for (kk = 0; kk < ILP_nneigh; kk++) {
      k = ILP_neigh[kk];
      if (k == i) continue;
      // derivatives of the product of rij and ni respect to rk, k=0,1,2, where atom k is the neighbors of atom i
      fk[0] = dnormdxk[0][0][kk] * dproddni[0] + dnormdxk[1][0][kk] * dproddni[1] +
          dnormdxk[2][0][kk] * dproddni[2];
      fk[1] = dnormdxk[0][1][kk] * dproddni[0] + dnormdxk[1][1][kk] * dproddni[1] +
          dnormdxk[2][1][kk] * dproddni[2];
      fk[2] = dnormdxk[0][2][kk] * dproddni[0] + dnormdxk[1][2][kk] * dproddni[1] +
          dnormdxk[2][2][kk] * dproddni[2];

      f[k][0] += fk[0];
      f[k][1] += fk[1];
      f[k][2] += fk[2];

      delki[0] = x[k][0] - x[i][0];
      delki[1] = x[k][1] - x[i][1];
      delki[2] = x[k][2] - x[i][2];

      if (VFLAG_EITHER) {
        ev_tally_xyz(k, i, nlocal, newton_pair, 0.0, 0.0, fk[0], fk[1], fk[2], delki[0], delki[1],
                     delki[2]);
      }
    }
    f[i][0] +=
        dnormdxi[0][0] * dproddni[0] + dnormdxi[1][0] * dproddni[1] + dnormdxi[2][0] * dproddni[2];
    f[i][1] +=
        dnormdxi[0][1] * dproddni[0] + dnormdxi[1][1] * dproddni[1] + dnormdxi[2][1] * dproddni[2];
    f[i][2] +=
        dnormdxi[0][2] * dproddni[0] + dnormdxi[1][2] * dproddni[1] + dnormdxi[2][2] * dproddni[2];
  }    // loop over ii
}

/* ----------------------------------------------------------------------
   Calculate the normals for one atom
------------------------------------------------------------------------- */
void PairILPGrapheneHBNOpt::calc_single_normal(int i, int *ILP_neigh, int nneigh, double *normal,
                                               double (*dnormdri)[3], double (*dnormdrk)[3][3])
{
  int cont, id, ip, m;
  double nn, delx, dely, delz, nn2;
  double pv12[3], pv31[3], pv23[3], n1[3], dni[3], dnn[3][3], vet[3][3], dpvdri[3][3];
  double dn1[3][3][3], dpv12[3][3][3], dpv23[3][3][3], dpv31[3][3][3];

  double **x = atom->x;

  for (id = 0; id < 3; id++) {
    pv12[id] = 0.0;
    pv31[id] = 0.0;
    pv23[id] = 0.0;
    n1[id] = 0.0;
    dni[id] = 0.0;
    normal[id] = 0.0;
    for (ip = 0; ip < 3; ip++) {
      vet[ip][id] = 0.0;
      dnn[ip][id] = 0.0;
      dpvdri[ip][id] = 0.0;
      dnormdri[ip][id] = 0.0;
      for (m = 0; m < 3; m++) {
        dpv12[ip][id][m] = 0.0;
        dpv31[ip][id][m] = 0.0;
        dpv23[ip][id][m] = 0.0;
        dn1[ip][id][m] = 0.0;
        dnormdrk[ip][id][m] = 0.0;
      }
    }
  }

  const double xtp = x[i][0];
  const double ytp = x[i][1];
  const double ztp = x[i][2];

  cont = 0;
  int j, *jlist = ILP_neigh;
  const int jnum = nneigh;
  for (int jj = 0; jj < jnum; jj++) {
    j = jlist[jj];
    j &= NEIGHMASK;

    delx = x[j][0] - xtp;
    dely = x[j][1] - ytp;
    delz = x[j][2] - ztp;
    vet[cont][0] = delx;
    vet[cont][1] = dely;
    vet[cont][2] = delz;
    cont++;
  }

  if (cont <= 1) {
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 1.0;
    for (id = 0; id < 3; id++) {
      for (ip = 0; ip < 3; ip++) {
        dnormdri[id][ip] = 0.0;
        for (m = 0; m < 3; m++) { dnormdrk[id][ip][m] = 0.0; }
      }
    }
  } else if (cont == 2) {
    pv12[0] = vet[0][1] * vet[1][2] - vet[1][1] * vet[0][2];
    pv12[1] = vet[0][2] * vet[1][0] - vet[1][2] * vet[0][0];
    pv12[2] = vet[0][0] * vet[1][1] - vet[1][0] * vet[0][1];

    // derivatives of pv12[0] to ri

    dpvdri[0][0] = 0.0;
    dpvdri[0][1] = vet[0][2] - vet[1][2];
    dpvdri[0][2] = vet[1][1] - vet[0][1];

    // derivatives of pv12[1] to ri

    dpvdri[1][0] = vet[1][2] - vet[0][2];
    dpvdri[1][1] = 0.0;
    dpvdri[1][2] = vet[0][0] - vet[1][0];

    // derivatives of pv12[2] to ri

    dpvdri[2][0] = vet[0][1] - vet[1][1];
    dpvdri[2][1] = vet[1][0] - vet[0][0];
    dpvdri[2][2] = 0.0;

    dpv12[0][0][0] = 0.0;
    dpv12[0][1][0] = vet[1][2];
    dpv12[0][2][0] = -vet[1][1];
    dpv12[1][0][0] = -vet[1][2];
    dpv12[1][1][0] = 0.0;
    dpv12[1][2][0] = vet[1][0];
    dpv12[2][0][0] = vet[1][1];
    dpv12[2][1][0] = -vet[1][0];
    dpv12[2][2][0] = 0.0;

    // derivatives respect to the second neighbor, atom l

    dpv12[0][0][1] = 0.0;
    dpv12[0][1][1] = -vet[0][2];
    dpv12[0][2][1] = vet[0][1];
    dpv12[1][0][1] = vet[0][2];
    dpv12[1][1][1] = 0.0;
    dpv12[1][2][1] = -vet[0][0];
    dpv12[2][0][1] = -vet[0][1];
    dpv12[2][1][1] = vet[0][0];
    dpv12[2][2][1] = 0.0;

    // derivatives respect to the third neighbor, atom n
    // derivatives of pv12 to rn is zero

    for (id = 0; id < 3; id++)
      for (ip = 0; ip < 3; ip++) dpv12[id][ip][2] = 0.0;

    n1[0] = pv12[0];
    n1[1] = pv12[1];
    n1[2] = pv12[2];

    // the magnitude of the normal vector
    nn2 = n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2];
    nn = sqrt(nn2);
    if (nn == 0) error->one(FLERR, "The magnitude of the normal vector is zero");

    // the unit normal vector

    normal[0] = n1[0] / nn;
    normal[1] = n1[1] / nn;
    normal[2] = n1[2] / nn;

    // derivatives of nn, dnn:3x1 vector

    dni[0] = (n1[0] * dpvdri[0][0] + n1[1] * dpvdri[1][0] + n1[2] * dpvdri[2][0]) / nn;
    dni[1] = (n1[0] * dpvdri[0][1] + n1[1] * dpvdri[1][1] + n1[2] * dpvdri[2][1]) / nn;
    dni[2] = (n1[0] * dpvdri[0][2] + n1[1] * dpvdri[1][2] + n1[2] * dpvdri[2][2]) / nn;

    // derivatives of unit vector ni respect to ri, the result is 3x3 matrix

    for (id = 0; id < 3; id++) {
      for (ip = 0; ip < 3; ip++) {
        dnormdri[id][ip] = dpvdri[id][ip] / nn - n1[id] * dni[ip] / nn2;
      }
    }

    // derivatives of non-normalized normal vector, dn1:3x3x3 array

    for (id = 0; id < 3; id++) {
      for (ip = 0; ip < 3; ip++) {
        for (m = 0; m < 3; m++) dn1[id][ip][m] = dpv12[id][ip][m];
      }
    }

    // derivatives of nn, dnn:3x3 vector
    // dnn[id][m]: the derivative of nn respect to r[id][m], id,m=0,1,2
    // r[id][m]: the id's component of atom m

    for (m = 0; m < 3; m++) {
      for (id = 0; id < 3; id++)
        dnn[id][m] = (n1[0] * dn1[0][id][m] + n1[1] * dn1[1][id][m] + n1[2] * dn1[2][id][m]) / nn;
    }

    // dnormdrk[id][ip][m][i]: the derivative of normal[id] respect to r[ip][m], id,ip=0,1,2
    // for atom m, which is a neighbor atom of atom i, m=0,jnum-1

    for (m = 0; m < 3; m++) {
      for (id = 0; id < 3; id++) {
        for (ip = 0; ip < 3; ip++)
          dnormdrk[id][ip][m] = dn1[id][ip][m] / nn - n1[id] * dnn[ip][m] / nn2;
      }
    }
  }

  //##############################################################################################

  else if (cont == 3) {
    pv12[0] = vet[0][1] * vet[1][2] - vet[1][1] * vet[0][2];
    pv12[1] = vet[0][2] * vet[1][0] - vet[1][2] * vet[0][0];
    pv12[2] = vet[0][0] * vet[1][1] - vet[1][0] * vet[0][1];

    // derivatives respect to the first neighbor, atom k

    dpv12[0][0][0] = 0.0;
    dpv12[0][1][0] = vet[1][2];
    dpv12[0][2][0] = -vet[1][1];
    dpv12[1][0][0] = -vet[1][2];
    dpv12[1][1][0] = 0.0;
    dpv12[1][2][0] = vet[1][0];
    dpv12[2][0][0] = vet[1][1];
    dpv12[2][1][0] = -vet[1][0];
    dpv12[2][2][0] = 0.0;

    // derivatives respect to the second neighbor, atom l

    dpv12[0][0][1] = 0.0;
    dpv12[0][1][1] = -vet[0][2];
    dpv12[0][2][1] = vet[0][1];
    dpv12[1][0][1] = vet[0][2];
    dpv12[1][1][1] = 0.0;
    dpv12[1][2][1] = -vet[0][0];
    dpv12[2][0][1] = -vet[0][1];
    dpv12[2][1][1] = vet[0][0];
    dpv12[2][2][1] = 0.0;

    // derivatives respect to the third neighbor, atom n

    for (id = 0; id < 3; id++) {
      for (ip = 0; ip < 3; ip++) dpv12[id][ip][2] = 0.0;
    }

    pv31[0] = vet[2][1] * vet[0][2] - vet[0][1] * vet[2][2];
    pv31[1] = vet[2][2] * vet[0][0] - vet[0][2] * vet[2][0];
    pv31[2] = vet[2][0] * vet[0][1] - vet[0][0] * vet[2][1];

    // derivatives respect to the first neighbor, atom k

    dpv31[0][0][0] = 0.0;
    dpv31[0][1][0] = -vet[2][2];
    dpv31[0][2][0] = vet[2][1];
    dpv31[1][0][0] = vet[2][2];
    dpv31[1][1][0] = 0.0;
    dpv31[1][2][0] = -vet[2][0];
    dpv31[2][0][0] = -vet[2][1];
    dpv31[2][1][0] = vet[2][0];
    dpv31[2][2][0] = 0.0;

    // derivatives respect to the third neighbor, atom n

    dpv31[0][0][2] = 0.0;
    dpv31[0][1][2] = vet[0][2];
    dpv31[0][2][2] = -vet[0][1];
    dpv31[1][0][2] = -vet[0][2];
    dpv31[1][1][2] = 0.0;
    dpv31[1][2][2] = vet[0][0];
    dpv31[2][0][2] = vet[0][1];
    dpv31[2][1][2] = -vet[0][0];
    dpv31[2][2][2] = 0.0;

    // derivatives respect to the second neighbor, atom l

    for (id = 0; id < 3; id++) {
      for (ip = 0; ip < 3; ip++) dpv31[id][ip][1] = 0.0;
    }

    pv23[0] = vet[1][1] * vet[2][2] - vet[2][1] * vet[1][2];
    pv23[1] = vet[1][2] * vet[2][0] - vet[2][2] * vet[1][0];
    pv23[2] = vet[1][0] * vet[2][1] - vet[2][0] * vet[1][1];

    // derivatives respect to the second neighbor, atom k

    for (id = 0; id < 3; id++) {
      for (ip = 0; ip < 3; ip++) { dpv23[id][ip][0] = 0.0; }
    }

    // derivatives respect to the second neighbor, atom l

    dpv23[0][0][1] = 0.0;
    dpv23[0][1][1] = vet[2][2];
    dpv23[0][2][1] = -vet[2][1];
    dpv23[1][0][1] = -vet[2][2];
    dpv23[1][1][1] = 0.0;
    dpv23[1][2][1] = vet[2][0];
    dpv23[2][0][1] = vet[2][1];
    dpv23[2][1][1] = -vet[2][0];
    dpv23[2][2][1] = 0.0;

    // derivatives respect to the third neighbor, atom n

    dpv23[0][0][2] = 0.0;
    dpv23[0][1][2] = -vet[1][2];
    dpv23[0][2][2] = vet[1][1];
    dpv23[1][0][2] = vet[1][2];
    dpv23[1][1][2] = 0.0;
    dpv23[1][2][2] = -vet[1][0];
    dpv23[2][0][2] = -vet[1][1];
    dpv23[2][1][2] = vet[1][0];
    dpv23[2][2][2] = 0.0;

    //############################################################################################
    // average the normal vectors by using the 3 neighboring planes

    n1[0] = (pv12[0] + pv31[0] + pv23[0]) / cont;
    n1[1] = (pv12[1] + pv31[1] + pv23[1]) / cont;
    n1[2] = (pv12[2] + pv31[2] + pv23[2]) / cont;

    // the magnitude of the normal vector

    nn2 = n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2];
    nn = sqrt(nn2);
    if (nn == 0) error->one(FLERR, "The magnitude of the normal vector is zero");

    // the unit normal vector

    normal[0] = n1[0] / nn;
    normal[1] = n1[1] / nn;
    normal[2] = n1[2] / nn;

    // for the central atoms, dnormdri is always zero

    for (id = 0; id < 3; id++) {
      for (ip = 0; ip < 3; ip++) dnormdri[id][ip] = 0.0;
    }

    // derivatives of non-normalized normal vector, dn1:3x3x3 array

    for (id = 0; id < 3; id++) {
      for (ip = 0; ip < 3; ip++) {
        for (m = 0; m < 3; m++)
          dn1[id][ip][m] = (dpv12[id][ip][m] + dpv23[id][ip][m] + dpv31[id][ip][m]) / cont;
      }
    }

    // derivatives of nn, dnn:3x3 vector
    // dnn[id][m]: the derivative of nn respect to r[id][m], id,m=0,1,2
    // r[id][m]: the id's component of atom m

    for (m = 0; m < 3; m++) {
      for (id = 0; id < 3; id++)
        dnn[id][m] = (n1[0] * dn1[0][id][m] + n1[1] * dn1[1][id][m] + n1[2] * dn1[2][id][m]) / nn;
    }

    // dnormdrk[id][ip][m][i]: the derivative of normal[id] respect to r[ip][m], id,ip=0,1,2
    // for atom m, which is a neighbor atom of atom i, m=0,jnum-1

    for (m = 0; m < 3; m++) {
      for (id = 0; id < 3; id++) {
        for (ip = 0; ip < 3; ip++)
          dnormdrk[id][ip][m] = dn1[id][ip][m] / nn - n1[id] * dnn[ip][m] / nn2;
      }
    }
  } else {
    error->one(FLERR, "There are too many neighbors for calculating normals");
  }
}

/* ------------------------------------------------------------------------ */

bool check_vdw(tagint itag, tagint jtag, double *xi, double *xj)
{
  if (itag > jtag) {
    if ((itag + jtag) % 2 == 0) return false;
  } else if (itag < jtag) {
    if ((itag + jtag) % 2 == 1) return false;
  } else {
    if (xj[2] < xi[2]) return false;
    if (xj[2] == xi[2] && xj[1] < xi[1]) return false;
    if (xj[2] == xi[2] && xj[1] == xi[1] && xj[0] < xi[0]) return false;
  }
  return true;
}

/* ------------------------------------------------------------------------ */

void PairILPGrapheneHBNOpt::update_internal_list()
{
  int jnum_sum = 0;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  tagint *tag = atom->tag;
  double **x = atom->x;
  for (int ii = 0; ii < inum; ii++) { jnum_sum += numneigh[ilist[ii]]; }
  if (inum > inum_max) {
    memory->destroy(num_intra);
    memory->destroy(num_inter);
    memory->destroy(num_vdw);
    memory->sfree(first_layered_neigh);
    //golden ratio grow
    inum_max = (int) ceil(inum / 0.618);
    memory->create(num_intra, inum_max, "PairILPGrapheneHBN:intra_layer_count");
    memory->create(num_inter, inum_max, "PairILPGrapheneHBN:inter_layer_count");
    memory->create(num_vdw, inum_max, "PairILPGrapheneHBN:vdw_count");
    first_layered_neigh = (int **) memory->smalloc(inum_max * sizeof(int *),
                                                   "PairILPGrapheneHBN:first_layered_neigh");
  }
  if (jnum_sum > jnum_max) {
    memory->destroy(layered_neigh);
    jnum_max = (int) ceil(jnum_sum / 0.618);
    memory->create(layered_neigh, jnum_max, "PairILPGrapheneHBN:layered_neigh");
  }

  double cut_intra = 0;
  for (int i = 0; i < nparams; i++)
    if (params[i].rcut > cut_intra) { cut_intra = params[i].rcut; }

  double cut_intra_listsq = (cut_intra + neighbor->skin) * (cut_intra + neighbor->skin);

  int total_neigh = 0;
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    tagint itag = tag[i];
    int jnum = numneigh[i];
    int *jlist = firstneigh[i];
    int *jlist_layered = first_layered_neigh[i] = layered_neigh + total_neigh;
    int ninter = 0, nintra = 0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      if (atom->molecule[j] == atom->molecule[i]) {
        double delx = x[i][0] - x[j][0];
        double dely = x[i][1] - x[j][1];
        double delz = x[i][2] - x[j][2];
        double rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < cut_intra_listsq) jlist_layered[nintra++] = j;
      }
    }
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      tagint jtag = tag[j];
      if (atom->molecule[j] != atom->molecule[i]) {
        if (check_vdw(itag, jtag, x[i], x[j])) jlist_layered[nintra + ninter++] = j;
      }
    }
    num_vdw[i] = ninter;
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      tagint jtag = tag[j];
      if (atom->molecule[j] != atom->molecule[i]) {
        if (!check_vdw(itag, jtag, x[i], x[j])) jlist_layered[nintra + ninter++] = j;
      }
    }
    num_intra[i] = nintra;
    num_inter[i] = ninter;
    total_neigh += nintra + ninter;
  }
}

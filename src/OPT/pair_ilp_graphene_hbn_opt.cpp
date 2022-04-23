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
   Contributing author: Wengen Ouyang (Tel Aviv University)
   e-mail: w.g.ouyang at gmail dot com

   This is a full version of the potential described in
   [Maaravi et al, J. Phys. Chem. C 121, 22826-22835 (2017)]
   The definition of normals are the same as that in
   [Kolmogorov & Crespi, Phys. Rev. B 71, 235415 (2005)]
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
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>
#include <map>

using namespace LAMMPS_NS;
using namespace InterLayer;

#define MAXLINE 1024
#define DELTA 4
#define PGDELTA 1
static const char cite_ilp[] =
    "ilp/graphene/hbn potential doi:10.1021/acs.nanolett.8b02848\n"
    "@Article{Ouyang2018\n"
    " author = {W. Ouyang, D. Mandelli, M. Urbakh, and O. Hod},\n"
    " title = {Nanoserpents: Graphene Nanoribbon Motion on Two-Dimensional Hexagonal Materials},\n"
    " journal = {Nano Letters},\n"
    " volume =  18,\n"
    " pages =   {6009}\n"
    " year =    2018,\n"
    "}\n\n";

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

/* ---------------------------------------------------------------------- */

PairILPGrapheneHBNOpt::PairILPGrapheneHBNOpt(LAMMPS *lmp) : PairILPGrapheneHBN(lmp)
{
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  if (lmp->citeme) {
    lmp->citeme->add(cite_ilp);
    lmp->citeme->add(cite_ilp_cur);
  }

  nextra = 2;
  pvector = new double[nextra];

  // initialize element to parameter maps
  params = nullptr;
  cutILPsq = nullptr;

  nmax = 0;
  maxlocal = 0;

  normal = nullptr;
  dnormal = nullptr;
  dnormdri = nullptr;

  // always compute energy offset
  offset_flag = 1;

  // turn on the taper function by default
  tap_flag = 1;
}

/* ---------------------------------------------------------------------- */

PairILPGrapheneHBNOpt::~PairILPGrapheneHBNOpt()
{
  memory->destroy(normal);
  memory->destroy(dnormal);
  memory->destroy(dnormdri);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(offset);
  }

  memory->destroy(elem2param);
  memory->destroy(cutILPsq);
  memory->sfree(params);
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

  // need a full neighbor list, including neighbors of ghosts
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
}

/* ---------------------------------------------------------------------- */
void PairILPGrapheneHBNOpt::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  //refactor and optimize for this pair style
  computeILP(eflag, vflag);
  return;
}

/* ----------------------------------------------------------------------
  refactor and optimize for this pair style
------------------------------------------------------------------------- */

void PairILPGrapheneHBNOpt::computeILP(int eflag, int vflag)
{
  pvector[0] = pvector[1] = 0.0;
  int i, j, ii, jj, inum, jnum, itype, jtype, k, l, kk, ll;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp, delx, dely, delz, erep, fpair;
  double rsq, r, Rcut, r2inv, r6inv, r8inv, Tap, dTap, Vilp, TSvdw, TSvdw2inv, fsum;
  int *ilist, *jlist, *numneigh, **firstneigh;

  erep = 0.0;
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  //Vars for calc_normal
  int id, ip, m;
  double nn, nn2;
  double pv12[3], pv31[3], pv23[3], n1[3], dni[3], dnn[3][3], vet[3][3], dpvdri[3][3];
  double dn1[3][3][3], dpv12[3][3][3], dpv23[3][3][3], dpv31[3][3][3];
  double normal[3], dnormal[3][3][3], dnormdri[3][3];
  //Vars for calc_Frep
  double prodnorm1, fkcx, fkcy, fkcz;
  double rhosq1, exp0, exp1;
  double frho1, Erep, rdsq1, fpair1;
  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double delkj[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};
  double fkk[3][3], ekk[3], vkk[3][6];

  //more for overflow
  int ilp_neigh[6];

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    int itypem = map[itype];
    itag = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int nilp = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      int jtypem = map[jtype];
      jtag = tag[j];

      // two-body interactions from full neighbor list, skip half of them

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq != 0 && rsq < cutILPsq[itypem][jtypem] && atom->molecule[i] == atom->molecule[j]) {
        ilp_neigh[nilp] = j;
        vet[nilp][0] = -delx;
        vet[nilp][1] = -dely;
        vet[nilp][2] = -delz;
        nilp++;
      }
    }
    //nilp;
    for (id = 0; id < 3; id++) {
      pv12[id] = 0.0;
      pv31[id] = 0.0;
      pv23[id] = 0.0;
      n1[id] = 0.0;
      dni[id] = 0.0;
      normal[id] = 0.0;
      for (ip = 0; ip < 3; ip++) {
        //vet[ip][id] = 0.0;
        dnn[ip][id] = 0.0;
        dpvdri[ip][id] = 0.0;
        dnormdri[ip][id] = 0.0;
        for (m = 0; m < 3; m++) {
          dpv12[m][ip][id] = 0.0;
          dpv31[m][ip][id] = 0.0;
          dpv23[m][ip][id] = 0.0;
          dn1[m][ip][id] = 0.0;
          dnormal[m][ip][id] = 0.0;
        }
      }
    }

    if (nilp <= 1) {
      normal[0] = 0.0;
      normal[1] = 0.0;
      normal[2] = 1.0;
      for (id = 0; id < 3; id++) {
        for (ip = 0; ip < 3; ip++) {
          dnormdri[id][ip] = 0.0;
          for (m = 0; m < 3; m++) { dnormal[m][id][ip] = 0.0; }
        }
      }
    } else if (nilp == 2) {
      cross_deriv(pv12, dpv12, vet, 0, 1, 2);
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
      for (m = 0; m < 3; m++) {
        for (id = 0; id < 3; id++) {
          for (ip = 0; ip < 3; ip++) { dn1[m][id][ip] = dpv12[m][id][ip]; }
        }
      }
      calc_dnormal(dnormal, dn1, n1, nn, nn2);

    } else if (nilp == 3) {
      cross_deriv(pv12, dpv12, vet, 0, 1, 2);
      cross_deriv(pv31, dpv31, vet, 2, 0, 1);
      cross_deriv(pv23, dpv23, vet, 1, 2, 0);

      n1[0] = (pv12[0] + pv31[0] + pv23[0]) / jnum;
      n1[1] = (pv12[1] + pv31[1] + pv23[1]) / jnum;
      n1[2] = (pv12[2] + pv31[2] + pv23[2]) / jnum;
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
        for (ip = 0; ip < 3; ip++) { dnormdri[id][ip] = 0.0; }
      }

      // derivatives of non-normalized normal vector, dn1:3x3x3 array
      for (m = 0; m < 3; m++) {
        for (id = 0; id < 3; id++) {
          for (ip = 0; ip < 3; ip++) {
            dn1[m][id][ip] = (dpv12[m][id][ip] + dpv23[m][id][ip] + dpv31[m][id][ip]) / jnum;
          }
        }
      }
      calc_dnormal(dnormal, dn1, n1, nn, nn2);
    } else {
      error->one(FLERR, "There are too many neighbors for calculating normals");
    }
    for (kk = 0; kk < nilp; kk++) {
      fkk[kk][0] = 0.0;
      fkk[kk][1] = 0.0;
      fkk[kk][2] = 0.0;
      vkk[kk][0] = 0.0;
      vkk[kk][1] = 0.0;
      vkk[kk][2] = 0.0;
      vkk[kk][3] = 0.0;
      vkk[kk][4] = 0.0;
      vkk[kk][5] = 0.0;
      ekk[kk] = 0.0;
    }
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = type[j];
      jtag = tag[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      // only include the interation between different layers
      if (rsq < cutsq[itype][jtype] && atom->molecule[i] != atom->molecule[j]) {

        int iparam_ij = elem2param[map[itype]][map[jtype]];

        //Param& p = params[iparam_ij];
        Param *p = &params[iparam_ij];

        r = sqrt(rsq);
        double rinv = 1.0 / r;

        // turn on/off taper function
        if (tap_flag) {
          Rcut = sqrt(cutsq[itype][jtype]);
          Tap = calc_Tap(r, Rcut);
          dTap = calc_dTap(r, Rcut);
        } else {
          Tap = 1.0;
          dTap = 0.0;
        }

        int vdwflag = 1;
        if (itag > jtag) {
          if ((itag + jtag) % 2 == 0) vdwflag = 0;
        } else if (itag < jtag) {
          if ((itag + jtag) % 2 == 1) vdwflag = 0;
        } else {
          if (x[j][2] < ztmp) vdwflag = 0;
          if (x[j][2] == ztmp && x[j][1] < ytmp) vdwflag = 0;
          if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) vdwflag = 0;
        }

        double fvdw = 0, evdw = 0;

        if (vdwflag) {
          r2inv = rinv * rinv;
          r6inv = r2inv * r2inv * r2inv;
          r8inv = r6inv * r2inv;

          TSvdw = 1.0 + exp(-p->d * (r / p->seff - 1.0));
          TSvdw2inv = 1 / (TSvdw * TSvdw);    //pow(TSvdw,-2.0);
          double vvdw = -p->C6 * r6inv / TSvdw;

          // derivatives
          fpair = -6.0 * p->C6 * r8inv / TSvdw +
              p->C6 * p->d / p->seff * (TSvdw - 1.0) * TSvdw2inv * r8inv * r;
          fsum = fpair * Tap - vvdw * dTap * rinv;

          fvdw = fsum;
          evdw = vvdw * Tap;
          pvector[0] += evdw;
        }

        // Calculate the transverse distance
        prodnorm1 = normal[0] * delx + normal[1] * dely + normal[2] * delz;
        rhosq1 = rsq - prodnorm1 * prodnorm1;    // rho_ij
        rdsq1 = rhosq1 * p->delta2inv;           // (rho_ij/delta)^2

        // store exponents
        exp0 = exp(-p->lambda * (r - p->z0));
        exp1 = exp(-rdsq1);

        frho1 = exp1 * p->C;
        Erep = 0.5 * p->epsilon + frho1;
        Vilp = exp0 * Erep;

        // derivatives
        fpair = p->lambda * exp0 * rinv * Erep;
        fpair1 = 2.0 * exp0 * frho1 * p->delta2inv;
        fsum = fpair + fpair1;
        // derivatives of the product of rij and ni, the result is a vector
        dprodnorm1[0] = dnormdri[0][0] * delx + dnormdri[1][0] * dely + dnormdri[2][0] * delz;
        dprodnorm1[1] = dnormdri[0][1] * delx + dnormdri[1][1] * dely + dnormdri[2][1] * delz;
        dprodnorm1[2] = dnormdri[0][2] * delx + dnormdri[1][2] * dely + dnormdri[2][2] * delz;
        fp1[0] = prodnorm1 * normal[0] * fpair1;
        fp1[1] = prodnorm1 * normal[1] * fpair1;
        fp1[2] = prodnorm1 * normal[2] * fpair1;
        fprod1[0] = prodnorm1 * dprodnorm1[0] * fpair1;
        fprod1[1] = prodnorm1 * dprodnorm1[1] * fpair1;
        fprod1[2] = prodnorm1 * dprodnorm1[2] * fpair1;

        fkcx = (delx * fsum - fp1[0]) * Tap - Vilp * dTap * delx * rinv;
        fkcy = (dely * fsum - fp1[1]) * Tap - Vilp * dTap * dely * rinv;
        fkcz = (delz * fsum - fp1[2]) * Tap - Vilp * dTap * delz * rinv;

        //This should be no use because fkcx need a lot of variables
        //fi + fj + sum(fk) = 0
        //-sum(fk) = fi + fj = -fprod*Tap
        //sum(fk) = fprod * Tap
        //fj = -fi - fprod*Tap
        double ftotx = fvdw * delx + fkcx;
        double ftoty = fvdw * dely + fkcy;
        double ftotz = fvdw * delz + fkcz;
        f[i][0] += ftotx - fprod1[0] * Tap;
        f[i][1] += ftoty - fprod1[1] * Tap;
        f[i][2] += ftotz - fprod1[2] * Tap;
        f[j][0] -= ftotx;
        f[j][1] -= ftoty;
        f[j][2] -= ftotz;

        // calculate the forces acted on the neighbors of atom i from atom j
        //ILP_neighs_i = ilp_neigh;
        for (kk = 0; kk < nilp; kk++) {
          k = ilp_neigh[kk];
          //if (k == i) continue;
          // derivatives of the product of rij and ni respect to rk, k=0,1,2, where atom k is the neighbors of atom i
          dprodnorm1[0] =
              dnormal[kk][0][0] * delx + dnormal[kk][1][0] * dely + dnormal[kk][2][0] * delz;
          dprodnorm1[1] =
              dnormal[kk][0][1] * delx + dnormal[kk][1][1] * dely + dnormal[kk][2][1] * delz;
          dprodnorm1[2] =
              dnormal[kk][0][2] * delx + dnormal[kk][1][2] * dely + dnormal[kk][2][2] * delz;
          fk[0] = (-prodnorm1 * dprodnorm1[0] * fpair1) * Tap;
          fk[1] = (-prodnorm1 * dprodnorm1[1] * fpair1) * Tap;
          fk[2] = (-prodnorm1 * dprodnorm1[2] * fpair1) * Tap;
          fkk[kk][0] += fk[0];
          fkk[kk][1] += fk[1];
          fkk[kk][2] += fk[2];
          delkj[0] = x[i][0] + vet[kk][0] - x[j][0];
          delkj[1] = x[i][1] + vet[kk][1] - x[j][1];
          delkj[2] = x[i][2] + vet[kk][2] - x[j][2];
          //if (evflag) ev_tally_xyz(k,j,nlocal,newton_pair,0.0,0.0,fk[0],fk[1],fk[2],delkj[0],delkj[1],delkj[2]);
          if (vflag_atom) {
            if (evflag)
              ev_tally_buffer(k, j, ekk + kk, vkk[kk], eatom + j, vatom[j], nlocal, newton_pair,
                              0.0, 0.0, fk[0], fk[1], fk[2], delkj[0], delkj[1], delkj[2]);
          } else {
            if (evflag)
              ev_tally_buffer(k, j, ekk + kk, vkk[kk], eatom + j, NULL, nlocal, newton_pair, 0.0,
                              0.0, fk[0], fk[1], fk[2], delkj[0], delkj[1], delkj[2]);
          }
        }
        erep = Tap * Vilp;
        if (eflag) pvector[1] += erep;    // = Tap*Vilp;
        //if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,fkcx,fkcy,fkcz,delx,dely,delz);
        if (evflag)
          ev_tally_xyz(i, j, nlocal, newton_pair, erep + evdw, 0.0, ftotx, ftoty, ftotz, delx, dely,
                       delz);
      }
    }    // loop over jj
    for (kk = 0; kk < nilp; kk++) {
      int k = ilp_neigh[kk];
      f[k][0] += fkk[kk][0];
      f[k][1] += fkk[kk][1];
      f[k][2] += fkk[kk][2];
      if (eflag_atom) eatom[k] += ekk[kk];
      if (vflag_atom) {
        vatom[k][0] += vkk[kk][0];
        vatom[k][1] += vkk[kk][1];
        vatom[k][2] += vkk[kk][2];
        vatom[k][3] += vkk[kk][3];
        vatom[k][4] += vkk[kk][4];
        vatom[k][5] += vkk[kk][5];
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

inline void PairILPGrapheneHBNOpt::cross_deriv(double *pv, double (*dpv)[3][3], double (*vet)[3],
                                               int j, int k, int l)
{
  pv[0] = vet[j][1] * vet[k][2] - vet[k][1] * vet[j][2];
  pv[1] = vet[j][2] * vet[k][0] - vet[k][2] * vet[j][0];
  pv[2] = vet[j][0] * vet[k][1] - vet[k][0] * vet[j][1];

  dpv[j][0][0] = 0.0;
  dpv[j][0][1] = vet[k][2];
  dpv[j][0][2] = -vet[k][1];
  dpv[j][1][0] = -vet[k][2];
  dpv[j][1][1] = 0.0;
  dpv[j][1][2] = vet[k][0];
  dpv[j][2][0] = vet[k][1];
  dpv[j][2][1] = -vet[k][0];
  dpv[j][2][2] = 0.0;

  dpv[k][0][0] = 0.0;
  dpv[k][0][1] = -vet[j][2];
  dpv[k][0][2] = vet[j][1];
  dpv[k][1][0] = vet[j][2];
  dpv[k][1][1] = 0.0;
  dpv[k][1][2] = -vet[j][0];
  dpv[k][2][0] = -vet[j][1];
  dpv[k][2][1] = vet[j][0];
  dpv[k][2][2] = 0.0;

  dpv[l][0][0] = 0.0;
  dpv[l][0][1] = 0.0;
  dpv[l][0][2] = 0.0;
  dpv[l][1][0] = 0.0;
  dpv[l][1][1] = 0.0;
  dpv[l][1][2] = 0.0;
  dpv[l][2][0] = 0.0;
  dpv[l][2][1] = 0.0;
  dpv[l][2][2] = 0.0;
}

inline void PairILPGrapheneHBNOpt::calc_dnormal(double (*dnormal)[3][3], double (*dn1)[3][3],
                                                double *n1, double nn, double nn2)
{
  double dnn[3][3];
  int m, id, ip;
  // derivatives of nn, dnn:3x3 vector
  // dnn[id][m]: the derivative of nn respect to r[id][m], id,m=0,1,2
  // r[id][m]: the id's component of atom m
  for (m = 0; m < 3; m++) {
    for (id = 0; id < 3; id++) {
      dnn[id][m] = (n1[0] * dn1[m][0][id] + n1[1] * dn1[m][1][id] + n1[2] * dn1[m][2][id]) / nn;
    }
  }
  // dnormal[m][id][ip]: the derivative of normal[id] respect to r[ip][m], id,ip=0,1,2
  // for atom m, which is a neighbor atom of atom i, m=0,jnum-1
  for (m = 0; m < 3; m++) {
    for (id = 0; id < 3; id++) {
      for (ip = 0; ip < 3; ip++) {
        dnormal[m][id][ip] = dn1[m][id][ip] / nn - n1[id] * dnn[ip][m] / nn2;
      }
    }
  }
}

inline void PairILPGrapheneHBNOpt::ev_tally_buffer(int i, int j, double *ei, double *vi, double *ej,
                                                   double *vj, int nlocal, int newton_pair,
                                                   double evdwl, double ecoul, double fx, double fy,
                                                   double fz, double delx, double dely, double delz)
{
  double evdwlhalf, ecoulhalf, epairhalf, v[6];
  if (eflag_either) {
    if (eflag_global) {
      if (newton_pair) {
        eng_vdwl += evdwl;
        eng_coul += ecoul;
      } else {
        evdwlhalf = 0.5 * evdwl;
        ecoulhalf = 0.5 * ecoul;
        if (i < nlocal) {
          eng_vdwl += evdwlhalf;
          eng_coul += ecoulhalf;
        }
        if (j < nlocal) {
          eng_vdwl += evdwlhalf;
          eng_coul += ecoulhalf;
        }
      }
    }
    if (eflag_atom) {
      epairhalf = 0.5 * (evdwl + ecoul);
      if (newton_pair || i < nlocal) *ei += epairhalf;
      if (newton_pair || j < nlocal) *ej += epairhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx * fx;
    v[1] = dely * fy;
    v[2] = delz * fz;
    v[3] = delx * fy;
    v[4] = delx * fz;
    v[5] = dely * fz;

    if (vflag_global) {
      if (newton_pair) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i < nlocal) {
          virial[0] += 0.5 * v[0];
          virial[1] += 0.5 * v[1];
          virial[2] += 0.5 * v[2];
          virial[3] += 0.5 * v[3];
          virial[4] += 0.5 * v[4];
          virial[5] += 0.5 * v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5 * v[0];
          virial[1] += 0.5 * v[1];
          virial[2] += 0.5 * v[2];
          virial[3] += 0.5 * v[3];
          virial[4] += 0.5 * v[4];
          virial[5] += 0.5 * v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_pair || i < nlocal) {
        vi[0] += 0.5 * v[0];
        vi[1] += 0.5 * v[1];
        vi[2] += 0.5 * v[2];
        vi[3] += 0.5 * v[3];
        vi[4] += 0.5 * v[4];
        vi[5] += 0.5 * v[5];
      }
      if (newton_pair || j < nlocal) {
        vj[0] += 0.5 * v[0];
        vj[1] += 0.5 * v[1];
        vj[2] += 0.5 * v[2];
        vj[3] += 0.5 * v[3];
        vj[4] += 0.5 * v[4];
        vj[5] += 0.5 * v[5];
      }
    }
  }
}

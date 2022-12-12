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
   Contributing authors: D.K. Ward (donward@sandia.gov) and X.W. Zhou (Sandia)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   The formulation for this work follows (a) D.G. Pettifor, et al., Mat.
   Sci. and Eng. A365, 2-13, (2004);(b) D.A. Murdick, et al., Phys.
   Rev. B 73, 045206 (2006);(c) D.G. Pettifor and I.I. Oleinik., Phys
   Rev. Lett. 84, 4124 (2000); (d) D.K. Ward, et al., Phys. Rev. B 85,
   115206 (2012).

   Copyright (2012) Sandia Corporation.  Under the terms of Contract DE-
   AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   rights in this software.

   pairbop v 1.0 comes with no warranty of any kind.  pairbop v 1.0 is a
   copyrighted code that is distributed free-of-charge, under the terms
   of the GNU Public License (GPL).  See "Open-Source Rules"
   https://www.lammps.org/open_source.html
------------------------------------------------------------------------- */

// uncomment define to enable writing table files for debugging
// #define LMP_BOP_WRITE_TABLES 1

#include "pair_bop.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "tabular_function.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathSpecial::square;
using MathSpecial::cube;


/* ---------------------------------------------------------------------- */

PairBOP::PairParameters::PairParameters()
{
  cutB = 0.0;
  cutBsq = 0.0;
  cutL = 0.0;
  cutLsq = 0.0;
  betaS = nullptr;
  betaP = nullptr;
  rep = nullptr;
  cphi = nullptr;
  bo = nullptr;
}
PairBOP::PairParameters::~PairParameters()
{
  delete betaS;
  delete betaP;
  delete rep;
  delete cphi;
  delete bo;
}

/* ---------------------------------------------------------------------- */

PairBOP::PairBOP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  ghostneigh = 1;
  allocated = 0;

  pairParameters = nullptr;
  tripletParameters = nullptr;
  bop_elements = nullptr;
  bop_masses = nullptr;
  bop_types = 0;

  pairlist1 = nullptr;
  pairlist2 = nullptr;
  triplelist = nullptr;
  neineilimit = -1;

  BOP_index = nullptr;
  BOP_index2 = nullptr;
  BOP_total = nullptr;
  BOP_total2 = nullptr;
  neigh_index = nullptr;
  neigh_index2 = nullptr;
  cos_index = nullptr;
  atomlimit = -1;
  neighlimit = -1;
  neighlimit2 = -1;

  map = nullptr;
  bytes = 0.0;
  pi_a = nullptr;
  pro_delta = nullptr;
  pi_delta = nullptr;
  pi_p = nullptr;
  pi_c = nullptr;
  pro = nullptr;
  sigma_delta = nullptr;
  sigma_c = nullptr;
  sigma_a = nullptr;
  sigma_f = nullptr;
  sigma_k = nullptr;
  small3 = nullptr;
  setflag = nullptr;
  cutsq = nullptr;
  cutghost = nullptr;

  bt_sg = nullptr;
  bt_pi = nullptr;
  sglimit = -1;
  pilimit = -1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairBOP::~PairBOP()
{
  if (allocated) {
    delete[] pairParameters;
    delete[] tripletParameters;
    memory->destroy(elem2param);
    memory->destroy(elem3param);
    memory->destroy(pi_a);
    memory->destroy(pro_delta);
    memory->destroy(pi_delta);
    memory->destroy(pi_p);
    memory->destroy(pi_c);
    memory->destroy(pro);
    memory->destroy(sigma_delta);
    memory->destroy(sigma_c);
    memory->destroy(sigma_a);
    memory->destroy(sigma_f);
    memory->destroy(sigma_k);
    memory->destroy(small3);
  }

  memory->destroy(setflag);
  memory->destroy(cutsq);
  memory->destroy(cutghost);

  memory->destroy(pairlist1);
  memory->destroy(pairlist2);
  memory->destroy(triplelist);
  memory->destroy(BOP_index);
  memory->destroy(BOP_total);
  memory->destroy(BOP_index2);
  memory->destroy(BOP_total2);
  memory->destroy(neigh_index);
  memory->destroy(neigh_index2);
  memory->destroy(cos_index);
  memory->destroy(bt_sg);
  memory->destroy(bt_pi);
  bytes = 0.0;

  if (bop_elements)
    for (int i = 0; i < bop_types; i++) delete[] bop_elements[i];
  delete[] bop_elements;
  delete[] bop_masses;
}

/* ---------------------------------------------------------------------- */

void PairBOP::compute(int eflag, int vflag)
{
  int i, ii, j, jj;
  int nlisti, *ilist;
  tagint i_tag,j_tag, itype, jtype;
  int temp_ij;
  double sigB_0, piB_0;
  double dpr1, dpr2, ftmp1, ftmp2, ftmp3, dE;

  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int *iilist = list->ilist;
  int **firstneigh = list->firstneigh;

  ev_init(eflag,vflag);

  gneigh();

  // For non on the fly calculations cos and derivatives
  // are calculated in advance and stored

  for (ii = 0; ii < nlocal; ii++) {
    i = iilist[ii];
    i_tag = tag[i];
    itype = map[type[i]];
    ilist = firstneigh[i];
    nlisti = BOP_total[i];
    for (jj = 0; jj < nlisti; jj++) {
      temp_ij = BOP_index[i] + jj;
      j = ilist[neigh_index[temp_ij]];
      j_tag = tag[j];
      if (i_tag > j_tag) {
        if ((i_tag+j_tag) % 2 == 0) continue;
      } else if (i_tag < j_tag) {
        if ((i_tag+j_tag) % 2 == 1) continue;
      } else {
        if (x[j][2] < x[i][2]) continue;
        if (x[j][2] == x[i][2] && x[j][1] < x[i][1]) continue;
        if (x[j][2] == x[i][2] && x[j][1] == x[i][1] && x[j][0] < x[i][0]) continue;
      }
      jtype = map[type[j]];
      int param_ij = elem2param[itype][jtype];
      sigB_0 = SigmaBo(ii,jj);
      if (pi_a[param_ij] == 0) {
        piB_0 = 0.0;
      } else {
        piB_0 = PiBo(ii,jj);
      }
      PairList1 & pl_ij = pairlist1[temp_ij];
      dpr1 = (pl_ij.dRep - 2.0*pl_ij.dBetaS*sigB_0 -
              2.0*pl_ij.dBetaP*piB_0) / pl_ij.r;
      ftmp1 = dpr1 * pl_ij.dis[0];
      ftmp2 = dpr1 * pl_ij.dis[1];
      ftmp3 = dpr1 * pl_ij.dis[2];
      f[i][0] += ftmp1;
      f[i][1] += ftmp2;
      f[i][2] += ftmp3;
      f[j][0] -= ftmp1;
      f[j][1] -= ftmp2;
      f[j][2] -= ftmp3;
      dE = pl_ij.rep - 2.0*pl_ij.betaS*sigB_0 - 2.0*pl_ij.betaP*piB_0;
      if (evflag) ev_tally(i,j,nlocal,newton_pair, dE, 0.0, -dpr1,
                           pl_ij.dis[0],pl_ij.dis[1],pl_ij.dis[2]);
    }
    nlisti = BOP_total2[i];
    for (jj = 0; jj < nlisti; jj++) {
      temp_ij = BOP_index2[i] + jj;
      j = ilist[neigh_index2[temp_ij]];
      j_tag = tag[j];
      if (i_tag > j_tag) {
        if ((i_tag+j_tag) % 2 == 0) continue;
      } else if (i_tag < j_tag) {
        if ((i_tag+j_tag) % 2 == 1) continue;
      } else {
        if (x[j][2] < x[i][2]) continue;
        if (x[j][2] == x[i][2] && x[j][1] < x[i][1]) continue;
        if (x[j][2] == x[i][2] && x[j][1] == x[i][1] && x[j][0] < x[i][0]) continue;
      }
      PairList2 & p2_ij = pairlist2[temp_ij];
      dpr2 = -p2_ij.dRep / p2_ij.r;
      ftmp1 = dpr2 * p2_ij.dis[0];
      ftmp2 = dpr2 * p2_ij.dis[1];
      ftmp3 = dpr2 * p2_ij.dis[2];
      f[i][0] += ftmp1;
      f[i][1] += ftmp2;
      f[i][2] += ftmp3;
      f[j][0] -= ftmp1;
      f[j][1] -= ftmp2;
      f[j][2] -= ftmp3;
      dE = -p2_ij.rep;
      if (evflag) ev_tally(i,j,nlocal,newton_pair, dE, 0.0, -dpr2,
                           p2_ij.dis[0],p2_ij.dis[1],p2_ij.dis[2]);
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBOP::allocate()
{
  allocated = 1;

  delete[] pairParameters;
  delete[] tripletParameters;
  memory->destroy(elem2param);
  memory->destroy(elem3param);

  memory->destroy(pi_a);
  memory->destroy(pro_delta);
  memory->destroy(pi_delta);
  memory->destroy(pi_p);
  memory->destroy(pi_c);
  memory->destroy(pro);
  memory->destroy(sigma_delta);
  memory->destroy(sigma_c);
  memory->destroy(sigma_a);
  memory->destroy(sigma_f);
  memory->destroy(sigma_k);
  memory->destroy(small3);

  pairParameters = new PairParameters[npairs];
  tripletParameters = new TabularFunction[ntriples];
  memory->create(elem2param,bop_types,bop_types,"BOP:elem2param");
  memory->create(elem3param,bop_types,bop_types,bop_types,"BOP:elem3param");
  bytes += (double)npairs*sizeof(PairParameters) +
    (double)ntriples*sizeof(TabularFunction) + square(bop_types)*sizeof(int) +
    cube(bop_types)*sizeof(int);

  memory->create(pi_a,npairs,"BOP:pi_a");
  memory->create(pro_delta,bop_types,"BOP:pro_delta");
  memory->create(pi_delta,npairs,"BOP:pi_delta");
  memory->create(pi_p,bop_types,"BOP:pi_p");
  memory->create(pi_c,npairs,"BOP:pi_c");
  memory->create(pro,bop_types,"BOP:pro");
  memory->create(sigma_delta,npairs,"BOP:sigma_delta");
  memory->create(sigma_c,npairs,"BOP:sigma_c");
  memory->create(sigma_a,npairs,"BOP:sigma_a");
  memory->create(sigma_f,npairs,"BOP:sigma_f");
  memory->create(sigma_k,npairs,"BOP:sigma_k");
  memory->create(small3,npairs,"BOP:small3");
  bytes += (9*npairs + 3*bop_types) * sizeof(double);

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBOP::settings(int narg, char **arg)
{
  otfly = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"save") == 0) {
      otfly = 0;
      iarg++;
    } else error->all(FLERR,"Illegal pair_style command");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs(Updated: D.K. Ward 05/06/10)
------------------------------------------------------------------------- */

void PairBOP::coeff(int narg, char **arg)
{
  const int np1 = atom->ntypes+1;
  delete[] map;
  map = new int[np1];
  memory->destroy(setflag);
  memory->destroy(cutsq);
  memory->destroy(cutghost);
  memory->create(setflag,np1,np1,"BOP:setflag");
  memory->create(cutsq,np1,np1,"BOP:cutsq");
  memory->create(cutghost,np1,np1,"BOP:cutghost");
  bytes = np1*np1*(sizeof (int) + 2.0*sizeof (double));

  map_element2type(narg-3, arg+3);

  // read the potential file

  read_table(arg[2]);

  // replace element indices in map with BOP potential indices
  // and check for missing elements

  if (comm->me == 0) {
    for (int i = 1; i < np1; i++) {
      int j;
      if (map[i] >= 0) {
        for (j = 0; j < bop_types; j++) {
          if (strcmp(elements[map[i]], bop_elements[j]) == 0) {
            map[i] = j;
            atom->set_mass(FLERR, i, bop_masses[j]);
            break;
          }
        }
        if (j == bop_types)
          error->one(FLERR,"Element {} not found in bop potential file {}",
                     elements[map[i]], arg[2]);
      }
    }
  }
  MPI_Bcast(map,np1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBOP::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style BOP requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style BOP requires newton pair on");

  if (utils::strmatch(force->pair_style,"^hybrid"))
    error->all(FLERR,"Pair style BOP is not compatible with hybrid pair styles");

  if ((neighbor->style == Neighbor::MULTI) || (neighbor->style == Neighbor::MULTI_OLD))
    error->all(FLERR,"Pair style BOP is not compatible with multi-cutoff neighbor lists");

  if (comm->mode != Comm::SINGLE)
    error->all(FLERR,"Pair style BOP is not compatible with multi-cutoff communication");

  // check that user sets comm->cutghostuser to 3x the max BOP cutoff

  if (comm->cutghostuser-0.001 < 3.0*cutmax)
    error->all(FLERR,"Pair style bop requires setting a communication cutoff "
               "of at least {:.4}",3.0*cutmax);

  // need a full neighbor list and neighbors of ghosts

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
}

/* ---------------------------------------------------------------------- */

double PairBOP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  int itype = map[i];
  int jtype = map[j];
  int param_ij = elem2param[itype][jtype];
  PairParameters & p_ij = pairParameters[param_ij];
  cutghost[i][j] = p_ij.cutB;
  if (cutghost[i][j] < p_ij.cutL) {
    cutghost[i][j] = p_ij.cutL;
  }
  cutsq[i][j] = cutghost[i][j]*cutghost[i][j];
  cutghost[j][i] = cutghost[i][j];
  cutsq[j][i] = cutsq[i][j];
  return cutghost[i][j];
}

/* ----------------------------------------------------------------------
   create BOP neighbor list from main neighbor list
   BOP neighbor list stores neighbors of ghost atoms
   BOP requires neighbor's of k if k is a neighbor of
   j and j is a neighbor of i
------------------------------------------------------------------------- */

void PairBOP::gneigh()
{
  int i, ii, j, jj, k, kk, temp_ij, temp_ik, temp_jik;
  int itype, jtype, ktype;
  int max_check, max_check2, neigh_total, neigh_total2, cos_total;
  double dis_ij[3], rsq_ij, r_ij;
  int *ilist, nlisti;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double **x = atom->x;
  int *type = atom->type;
  int *iilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  if (BOP_index) {
    if (atomlimit < nall) {
      bytes += 4*(nall-atomlimit) * sizeof (int);
      if (!otfly) bytes += (nall-atomlimit) * sizeof (int);
      atomlimit = nall;
      memory->grow(BOP_index,atomlimit,"BOP:BOP_index");
      memory->grow(BOP_index2,atomlimit,"BOP:BOP_index2");
      memory->grow(BOP_total,atomlimit,"BOP:BOP_total");
      memory->grow(BOP_total2,atomlimit,"BOP:BOP_total2");
      if (!otfly) memory->grow(cos_index,atomlimit,"BOP:cos_index");
    }
  } else {
    atomlimit = nall;
    memory->create(BOP_index,atomlimit,"BOP:BOP_index");
    memory->create(BOP_index2,atomlimit,"BOP:BOP_index2");
    memory->create(BOP_total,atomlimit,"BOP:BOP_total");
    memory->create(BOP_total2,atomlimit,"BOP:BOP_total");
    bytes += 4*atomlimit * sizeof (int);
    if (!otfly) {
      memory->create(cos_index,atomlimit,"BOP:cos_index");
      bytes += atomlimit * sizeof (int);
    }
  }

  neigh_total=0;
  neigh_total2=0;
  for (ii = 0; ii < nall; ii++) {
    if (ii < nlocal)
      i = iilist[ii];
    else
      i=ii;
    BOP_index[i] = neigh_total;
    BOP_index2[i] = neigh_total2;
    itype = map[type[i]];
    ilist = firstneigh[i];
    max_check = 0 ;
    max_check2 = 0;
    for (jj = 0; jj < numneigh[i]; jj++) {
      j = ilist[jj];
      jtype = map[type[j]];
      int param_ij = elem2param[itype][jtype];
      PairParameters & p_ij = pairParameters[param_ij];
      dis_ij[0] = x[j][0] - x[i][0];
      dis_ij[1] = x[j][1] - x[i][1];
      dis_ij[2] = x[j][2] - x[i][2];
      rsq_ij = dis_ij[0]*dis_ij[0] + dis_ij[1]*dis_ij[1] + dis_ij[2]*dis_ij[2];
      if (rsq_ij < p_ij.cutBsq) {
        r_ij = sqrt(rsq_ij);
        temp_ij = neigh_total + max_check;
        if (pairlist1) {
          if (neighlimit <= temp_ij) {
            neighlimit += 2*nall;
            memory->grow(pairlist1,neighlimit,"BOP:pairlist1");
            memory->grow(neigh_index,neighlimit,"BOP:neigh_index");
            bytes += 2*nall * (sizeof(PairList1)+sizeof (int));
          }
        } else {
          neighlimit = 12*nall;
          memory->create(pairlist1,neighlimit,"BOP:pairlist1");
          memory->create(neigh_index,neighlimit,"BOP:neigh_index");
          bytes += 12*nall * (sizeof(PairList1)+sizeof (int));
        }
        neigh_index[temp_ij] = jj;
        PairList1 & pl_ij = pairlist1[temp_ij];
        pl_ij.r = r_ij;
        pl_ij.dis[0] = dis_ij[0];
        pl_ij.dis[1] = dis_ij[1];
        pl_ij.dis[2] = dis_ij[2];
        (p_ij.betaS)->value(r_ij, pl_ij.betaS, 1, pl_ij.dBetaS, 1);
        if (pi_a[param_ij] != 0) {
          (p_ij.betaP)->value(r_ij, pl_ij.betaP, 1, pl_ij.dBetaP, 1);
        }
        (p_ij.rep)->value(r_ij, pl_ij.rep, 1, pl_ij.dRep, 1);
        max_check++;
      }
      if (rsq_ij < p_ij.cutLsq) {
        r_ij = sqrt(rsq_ij);
        temp_ij = neigh_total2 + max_check2;
        if (pairlist2) {
          if (neighlimit2 <= temp_ij) {
            neighlimit2 += 2*nall;
            memory->grow(pairlist2,neighlimit2,"BOP:pairlist2");
            memory->grow(neigh_index2,neighlimit2,"BOP:neigh_index2");
            bytes += 2*nall * (sizeof(PairList2)+sizeof (int));
          }
        } else {
          neighlimit2 = 12*nall;
          memory->create(pairlist2,neighlimit2,"BOP:pairlist2");
          memory->create(neigh_index2,neighlimit2,"BOP:neigh_index2");
          bytes += 12*nall * (sizeof(PairList2)+sizeof (int));
        }
        neigh_index2[temp_ij] = jj;
        PairList2 & pl_ij = pairlist2[temp_ij];
        pl_ij.r = r_ij;
        pl_ij.dis[0] = dis_ij[0];
        pl_ij.dis[1] = dis_ij[1];
        pl_ij.dis[2] = dis_ij[2];
        (p_ij.cphi)->value(r_ij, pl_ij.rep, 1, pl_ij.dRep, 1);
        max_check2++;
      }
    }
    BOP_total[i] = max_check;
    BOP_total2[i] = max_check2;
    neigh_total += max_check;
    neigh_total2 += max_check2;
  }

  if (!otfly) {
    cos_total = 0;
    for (ii = 0; ii < nall; ii++) {
      if (ii < nlocal)
        i = iilist[ii];
      else
        i=ii;
      cos_index[i] = cos_total;
      itype = map[type[i]];
      nlisti = BOP_total[i];
      ilist = firstneigh[i];
      max_check = 0;
      for (jj = 0; jj < nlisti; jj++) {
        temp_ij = BOP_index[i] + jj;
        j = ilist[neigh_index[temp_ij]];
        jtype = map[type[j]];
        PairList1 & pl_ij = pairlist1[temp_ij];
        for (kk = jj + 1; kk < nlisti; kk++) {
          temp_ik = BOP_index[i] + kk;
          k = ilist[neigh_index[temp_ik]];
          ktype = map[type[k]];
          int param_jik = elem3param[jtype][itype][ktype];
          PairList1 & pl_ik = pairlist1[temp_ik];
          temp_jik = cos_index[i] + max_check;
          if (triplelist) {
            if (neineilimit <= temp_jik) {
              neineilimit += 6*nall;
              memory->grow(triplelist,neineilimit,"BOP:triplelist");
              bytes += 6*nall * sizeof(TripleList);
            }
          } else {
            neineilimit = 6*11*nall;
            memory->create(triplelist,neineilimit,"BOP:triplelist");
            bytes += 6*11*nall * sizeof(TripleList);
          }
          TripleList & tl_jik = triplelist[temp_jik];
          angle(pl_ij.r, pl_ij.dis, pl_ik.r, pl_ik.dis, tl_jik.cosAng,
                tl_jik.dCosAngj, tl_jik.dCosAngk);
          tripletParameters[param_jik].value(tl_jik.cosAng, tl_jik.G, 1, tl_jik.dG, 1);
          max_check++;
        }
      }
      cos_total += max_check;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairBOP::angle(double r1, double *dis1, double r2, double *dis2,
                    double &ang, double *dAngj, double *dAngk)
{
  double rj1k1, rj2k1, rj1k2;
  rj1k1 = r1 * r2;
  rj2k1 = rj1k1 * r1;
  rj1k2 = rj1k1 * r2;
  ang = (dis1[0]*dis2[0] + dis1[1]*dis2[1] + dis1[2]*dis2[2]) / rj1k1;
  dAngj[0] = (dis2[0]*r1 - ang*dis1[0]*r2) / rj2k1;
  dAngj[1] = (dis2[1]*r1 - ang*dis1[1]*r2) / rj2k1;
  dAngj[2] = (dis2[2]*r1 - ang*dis1[2]*r2) / rj2k1;
  dAngk[0] = (dis1[0]*r2 - ang*dis2[0]*r1) / rj1k2;
  dAngk[1] = (dis1[1]*r2 - ang*dis2[1]*r1) / rj1k2;
  dAngk[2] = (dis1[2]*r2 - ang*dis2[2]*r1) / rj1k2;
}

/* ---------------------------------------------------------------------- */

/*  The formulation differs slightly to avoid negative square roots
    in the calculation of Sigma^(1/2) of (a) Eq. 6 and (b) Eq. 11 */

/* ---------------------------------------------------------------------- */

/*  The formulation differs slightly to avoid negative square roots
    in the calculation of Theta_pi,ij of (a) Eq. 36 and (b) Eq. 18 */

double PairBOP::SigmaBo(int itmp, int jtmp)
{
  double sigB, sigB1, dsigB, ftmp[3], xtmp[3];
  int i, j, k;
  int itype, jtype, ktype;
  int nb_t, nb_ij, nb_ik, nb_jk, bt_i, bt_j;
  int n_ji, n_jk, n_ki, n_kj, n_jik, n_ijk, n_ikj, pass_jk;
  int temp_ij, temp_ik, temp_jk, temp_kk, temp_jik, temp_ijk, temp_ikj;
  int *ilist, *jlist, *klist;
  int nlisti, nlistj, nlistk;
  double r_ij, dis_ij[3], r_ik, dis_ik[3], r_ji, dis_ji[3],
    r_jk, dis_jk[3], r_ki, dis_ki[3], r_kj, dis_kj[3];
  double dCosAngj[3], dCosAngk[3];
  double *gj, *gk;
  double betaS_ij, dBetaS_ij, betaS_ik, dBetaS_ik, betaS_jk, dBetaS_jk;
  double cosAng_jik, dcA_jik[3][2], cosAng_ijk, dcA_ijk[3][2],
    cosAng_ikj, dcA_ikj[3][2];
  int nfound, loop, temp_loop, nei_loop, nei;
  double AA, BB, EE1, FF, AAC, pp1;
  double gfactor1, gprime1, gfactor2, gprime2, gfactor3, gprime3,
    gfactor, gfactorsq, gsqprime, gcm1, gcm2, gcm3,
    rfactor, rcm1, rcm2;
  double agpdpr1, app1, bndtmp, bndtmp0, bndtmp1, bndtmp2,
    part0, part1, part2, part3, part4;
  int ktmp;

  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int *iilist = list->ilist;
  int **firstneigh = list->firstneigh;

  sigB = 0.0;
  if (itmp < nlocal) {
    i = iilist[itmp];
  } else {
    i=itmp;
  }

  nb_t = 0;
  memory_sg(nb_t);
  initial_sg(nb_t);

  itype = map[type[i]];
  ilist = firstneigh[i];
  nlisti = BOP_total[i];
  temp_ij = BOP_index[i] + jtmp;
  PairList1 & pl_ij = pairlist1[temp_ij];
  j = ilist[neigh_index[temp_ij]];
  jtype = map[type[j]];
  jlist = firstneigh[j];
  nlistj = BOP_total[j];
  int param_ij = elem2param[itype][jtype];
  PairParameters & p_ij = pairParameters[param_ij];
  nb_ij = nb_t;
  bt_sg[nb_ij].i = i;
  bt_sg[nb_ij].j = j;
  bt_sg[nb_ij].temp = temp_ij;
  nb_t++;
  memory_sg(nb_t);
  initial_sg(nb_t);

  for (loop = 0; loop < nlistj; loop++) {
    temp_loop = BOP_index[j] + loop;
    nei_loop = neigh_index[temp_loop];
    nei = jlist[nei_loop];
    if (x[nei][0]==x[i][0] && x[nei][1]==x[i][1] && x[nei][2]==x[i][2]) {
      n_ji = loop;
      break;
    }
  }

  dis_ij[0] = pl_ij.dis[0];
  dis_ij[1] = pl_ij.dis[1];
  dis_ij[2] = pl_ij.dis[2];
  r_ij = pl_ij.r;
  betaS_ij = pl_ij.betaS;
  dBetaS_ij = pl_ij.dBetaS;
  dis_ji[0] = -dis_ij[0];
  dis_ji[1] = -dis_ij[1];
  dis_ji[2] = -dis_ij[2];
  r_ji = r_ij;

  // AA-EE1 are the components making up Eq. 30 (a)

  AA = 0.0;
  BB = 0.0;
  EE1 = 0.0;

  // FF is the Beta_sigma^2 term

  FF = betaS_ij * betaS_ij;
  if (FF <= 0.000001) return(sigB);

  // agpdpr1 is derivative of FF w.r.t. r_ij

  agpdpr1 = 2.0 * betaS_ij * dBetaS_ij / r_ij;

  // dXX derivatives are taken with respect to all pairs contributing to the energy
  // nb_ij is derivative w.r.t. ij pair

  bt_sg[nb_ij].dFF[0] = agpdpr1 * dis_ij[0];
  bt_sg[nb_ij].dFF[1] = agpdpr1 * dis_ij[1];
  bt_sg[nb_ij].dFF[2] = agpdpr1 * dis_ij[2];

  // k is loop over all neighbors of i again with j neighbor of i

  for (ktmp = 0; ktmp < nlisti; ktmp++) {
    if (ktmp == jtmp) continue;
    temp_ik = BOP_index[i] + ktmp;
    PairList1 & pl_ik = pairlist1[temp_ik];
    k = ilist[neigh_index[temp_ik]];
    ktype = map[type[k]];
    klist = firstneigh[k];
    nlistk = BOP_total[k];

    dis_ik[0] = pl_ik.dis[0];
    dis_ik[1] = pl_ik.dis[1];
    dis_ik[2] = pl_ik.dis[2];
    r_ik = pl_ik.r;
    betaS_ik = pl_ik.betaS;
    dBetaS_ik = pl_ik.dBetaS;
    dis_ki[0] = -dis_ik[0];
    dis_ki[1] = -dis_ik[1];
    dis_ki[2] = -dis_ik[2];
    r_ki = r_ik;

    // find neighbors of k that are equal to i or j

    nfound = 0;
    pass_jk = 0;
    for (loop = 0; loop < nlistk; loop++) {
      temp_loop = BOP_index[k] + loop;
      nei_loop = neigh_index[temp_loop];
      nei = klist[nei_loop];
      if (x[nei][0]==x[i][0] && x[nei][1]==x[i][1] && x[nei][2]==x[i][2]) {
        n_ki = loop;
        nfound++;
        if (nfound == 2) break;
      }
      if (x[nei][0]==x[j][0] && x[nei][1]==x[j][1] && x[nei][2]==x[j][2]) {
        n_kj = loop;
        pass_jk = 1;
        nfound++;
        if (nfound == 2) break;
      }
    }

    nb_ik = nb_t;
    bt_sg[nb_ik].i = i;
    bt_sg[nb_ik].j = k;
    bt_sg[nb_ik].temp = temp_ik;
    nb_t++;
    memory_sg(nb_t);
    initial_sg(nb_t);
    if (pass_jk) {
      for (loop = 0; loop < nlistj; loop++) {
        temp_loop = BOP_index[j] + loop;
        nei_loop = neigh_index[temp_loop];
        nei = jlist[nei_loop];
        if (x[nei][0]==x[k][0] && x[nei][1]==x[k][1] && x[nei][2]==x[k][2]) {
          temp_jk = temp_loop;
          n_jk = loop;
          break;
        }
      }
      nb_jk = nb_t;
      bt_sg[nb_jk].i = j;
      bt_sg[nb_jk].j = k;
      bt_sg[nb_jk].temp = temp_jk;
      nb_t++;
      memory_sg(nb_t);
      initial_sg(nb_t);
    }

    if (!otfly) {
      if (jtmp < ktmp) {
        n_jik = jtmp*(2*nlisti-jtmp-1)/2 + (ktmp-jtmp)-1;
      } else {
        n_jik = ktmp*(2*nlisti-ktmp-1)/2 + (jtmp-ktmp)-1;
      }
      temp_jik = cos_index[i] + n_jik;
      TripleList & tl_jik = triplelist[temp_jik];
      if (jtmp < ktmp) {
        gj = tl_jik.dCosAngj;
        gk = tl_jik.dCosAngk;
      } else {
        gj = tl_jik.dCosAngk;
        gk = tl_jik.dCosAngj;
      }
      cosAng_jik = tl_jik.cosAng;
      dcA_jik[0][0] = gj[0];
      dcA_jik[1][0] = gj[1];
      dcA_jik[2][0] = gj[2];
      dcA_jik[0][1] = gk[0];
      dcA_jik[1][1] = gk[1];
      dcA_jik[2][1] = gk[2];
      gfactor1 = tl_jik.G;
      gprime1 = tl_jik.dG;
    } else {
      angle(r_ij, dis_ij, r_ik, dis_ik, cosAng_jik,
            dCosAngj, dCosAngk);
      int param_jik = elem3param[jtype][itype][ktype];
      tripletParameters[param_jik].value(cosAng_jik, gfactor1, 1, gprime1, 1);
      dcA_jik[0][0] = dCosAngj[0];
      dcA_jik[1][0] = dCosAngj[1];
      dcA_jik[2][0] = dCosAngj[2];
      dcA_jik[0][1] = dCosAngk[0];
      dcA_jik[1][1] = dCosAngk[1];
      dcA_jik[2][1] = dCosAngk[2];
    }
    gfactorsq = gfactor1*gfactor1;
    gsqprime = 2.0*gfactor1*gprime1;

    // AA is Eq. 34 (a) or Eq. 10 (c) for the i atom
    // 1st CC is Eq. 11 (c) for i atom where j & k=neighbor of i

    AA += gfactorsq * betaS_ik * betaS_ik;

    // agpdpr1 is derivative of AA w.r.t. rik
    // app1 is derivative of AA w.r.t. cos(theta_jik)

    agpdpr1 = 2.0 * gfactorsq * betaS_ik * dBetaS_ik / r_ik;
    app1 = betaS_ik * betaS_ik * gsqprime;
    bt_sg[nb_ij].dAA[0] += app1*dcA_jik[0][0];
    bt_sg[nb_ij].dAA[1] += app1*dcA_jik[1][0];
    bt_sg[nb_ij].dAA[2] += app1*dcA_jik[2][0];
    bt_sg[nb_ik].dAA[0] += app1*dcA_jik[0][1] + agpdpr1*dis_ik[0];
    bt_sg[nb_ik].dAA[1] += app1*dcA_jik[1][1] + agpdpr1*dis_ik[1];
    bt_sg[nb_ik].dAA[2] += app1*dcA_jik[2][1] + agpdpr1*dis_ik[2];

    // k' is loop over neighbors all neighbors of j with k a neighbor
    // of i and j a neighbor of i and determine which k' is k

    if (sigma_f[param_ij] == 0.5 || !sigma_k[param_ij] || !pass_jk) continue;
    PairList1 & pl_jk = pairlist1[temp_jk];
    dis_jk[0] = pl_jk.dis[0];
    dis_jk[1] = pl_jk.dis[1];
    dis_jk[2] = pl_jk.dis[2];
    r_jk = pl_jk.r;
    betaS_jk = pl_jk.betaS;
    dBetaS_jk = pl_jk.dBetaS;
    dis_kj[0] = -dis_jk[0];
    dis_kj[1] = -dis_jk[1];
    dis_kj[2] = -dis_jk[2];
    r_kj = r_jk;

    if (!otfly) {
      if (n_ji < n_jk) {
        n_ijk = n_ji*(2*nlistj-n_ji-1)/2 + (n_jk-n_ji)-1;
      } else {
        n_ijk = n_jk*(2*nlistj-n_jk-1)/2 + (n_ji-n_jk)-1;
      }
      temp_ijk = cos_index[j] + n_ijk;
      TripleList & tl_ijk = triplelist[temp_ijk];
      if (n_ji < n_jk) {
        gj = tl_ijk.dCosAngj;
        gk = tl_ijk.dCosAngk;
      } else {
        gj = tl_ijk.dCosAngk;
        gk = tl_ijk.dCosAngj;
      }
      cosAng_ijk = tl_ijk.cosAng;
      dcA_ijk[0][0] = gj[0];
      dcA_ijk[1][0] = gj[1];
      dcA_ijk[2][0] = gj[2];
      dcA_ijk[0][1] = gk[0];
      dcA_ijk[1][1] = gk[1];
      dcA_ijk[2][1] = gk[2];
      gfactor2 = tl_ijk.G;
      gprime2 = tl_ijk.dG;
      if (n_ki < n_kj) {
        n_ikj = n_ki*(2*nlistk-n_ki-1)/2 + (n_kj-n_ki)-1;
      } else {
        n_ikj = n_kj*(2*nlistk-n_kj-1)/2 + (n_ki-n_kj)-1;
      }
      temp_ikj = cos_index[k] + n_ikj;
      TripleList & tl_ikj = triplelist[temp_ikj];
      if (n_ki < n_kj) {
        gj = tl_ikj.dCosAngj;
        gk = tl_ikj.dCosAngk;
      } else {
        gj = tl_ikj.dCosAngk;
        gk = tl_ikj.dCosAngj;
      }
      cosAng_ikj = tl_ikj.cosAng;
      dcA_ikj[0][0] = gj[0];
      dcA_ikj[1][0] = gj[1];
      dcA_ikj[2][0] = gj[2];
      dcA_ikj[0][1] = gk[0];
      dcA_ikj[1][1] = gk[1];
      dcA_ikj[2][1] = gk[2];
      gfactor3 = tl_ikj.G;
      gprime3 = tl_ikj.dG;
    } else {
      angle(r_ji, dis_ji, r_jk, dis_jk, cosAng_ijk,
            dCosAngj, dCosAngk);
      int param_ijk = elem3param[itype][jtype][ktype];
      tripletParameters[param_ijk].value(cosAng_ijk, gfactor2, 1, gprime2, 1);
      dcA_ijk[0][0] = dCosAngj[0];
      dcA_ijk[1][0] = dCosAngj[1];
      dcA_ijk[2][0] = dCosAngj[2];
      dcA_ijk[0][1] = dCosAngk[0];
      dcA_ijk[1][1] = dCosAngk[1];
      dcA_ijk[2][1] = dCosAngk[2];
      angle(r_ki, dis_ki, r_kj, dis_kj, cosAng_ikj,
            dCosAngj, dCosAngk);
      int param_ikj = elem3param[itype][ktype][jtype];
      tripletParameters[param_ikj].value(cosAng_ikj, gfactor3, 1, gprime3, 1);
      dcA_ikj[0][0] = dCosAngj[0];
      dcA_ikj[1][0] = dCosAngj[1];
      dcA_ikj[2][0] = dCosAngj[2];
      dcA_ikj[0][1] = dCosAngk[0];
      dcA_ikj[1][1] = dCosAngk[1];
      dcA_ikj[2][1] = dCosAngk[2];
    }
    gfactor = gfactor1 * gfactor2 * gfactor3;
    rfactor = betaS_ik * betaS_jk;

    // EE1 is (b) Eq. 12

    EE1 += gfactor*rfactor;

    // rcm1 is derivative of EE1 w.r.t r_ik
    // rcm2 is derivative of EE1 w.r.t r_jk
    // gcm1 is derivative of EE1 w.r.t cos(theta_jik)
    // gcm2 is derivative of EE1 w.r.t cos(theta_ijk)
    // gcm3 is derivative of EE1 w.r.t cos(theta_ikj)

    rcm1 = gfactor*betaS_jk*dBetaS_ik / r_ik;
    rcm2 = gfactor*betaS_ik*dBetaS_jk / r_jk;
    gcm1 = rfactor*gprime1*gfactor2*gfactor3;
    gcm2 = rfactor*gfactor1*gprime2*gfactor3;
    gcm3 = rfactor*gfactor1*gfactor2*gprime3;
    bt_sg[nb_ij].dEE1[0] += gcm1*dcA_jik[0][0] - gcm2*dcA_ijk[0][0];
    bt_sg[nb_ij].dEE1[1] += gcm1*dcA_jik[1][0] - gcm2*dcA_ijk[1][0];
    bt_sg[nb_ij].dEE1[2] += gcm1*dcA_jik[2][0] - gcm2*dcA_ijk[2][0];
    bt_sg[nb_ik].dEE1[0] += gcm1*dcA_jik[0][1] + rcm1*dis_ik[0] -
      gcm3*dcA_ikj[0][0];
    bt_sg[nb_ik].dEE1[1] += gcm1*dcA_jik[1][1] + rcm1*dis_ik[1] -
      gcm3*dcA_ikj[1][0];
    bt_sg[nb_ik].dEE1[2] += gcm1*dcA_jik[2][1] + rcm1*dis_ik[2] -
      gcm3*dcA_ikj[2][0];
    bt_sg[nb_jk].dEE1[0] += gcm2*dcA_ijk[0][1] + rcm2*dis_jk[0] -
      gcm3*dcA_ikj[0][1];
    bt_sg[nb_jk].dEE1[1] += gcm2*dcA_ijk[1][1] + rcm2*dis_jk[1] -
      gcm3*dcA_ikj[1][1];
    bt_sg[nb_jk].dEE1[2] += gcm2*dcA_ijk[2][1] + rcm2*dis_jk[2] -
      gcm3*dcA_ikj[2][1];
  }

  for (ktmp = 0; ktmp < nlistj; ktmp++) {
    if (ktmp == n_ji) continue;
    temp_jk = BOP_index[j] + ktmp;
    PairList1 & pl_jk = pairlist1[temp_jk];
    n_jk = ktmp;
    k = jlist[neigh_index[temp_jk]];
    ktype = map[type[k]];
    dis_jk[0] = pl_jk.dis[0];
    dis_jk[1] = pl_jk.dis[1];
    dis_jk[2] = pl_jk.dis[2];
    r_jk = pl_jk.r;
    betaS_jk = pl_jk.betaS;
    dBetaS_jk = pl_jk.dBetaS;

    nb_jk = nb_t;
    bt_sg[nb_jk].i = j;
    bt_sg[nb_jk].j = k;
    bt_sg[nb_jk].temp = temp_jk;
    nb_t++;
    memory_sg(nb_t);
    initial_sg(nb_t);

    if (!otfly) {
      if (n_ji < n_jk) {
        n_ijk = n_ji*(2*nlistj-n_ji-1)/2 + (n_jk-n_ji)-1;
      } else {
        n_ijk = n_jk*(2*nlistj-n_jk-1)/2 + (n_ji-n_jk)-1;
      }
      temp_ijk = cos_index[j] + n_ijk;
      TripleList & tl_ijk = triplelist[temp_ijk];
      if (n_ji < n_jk) {
        gj = tl_ijk.dCosAngj;
        gk = tl_ijk.dCosAngk;
      } else {
        gj = tl_ijk.dCosAngk;
        gk = tl_ijk.dCosAngj;
      }
      cosAng_ijk = tl_ijk.cosAng;
      dcA_ijk[0][0] = gj[0];
      dcA_ijk[1][0] = gj[1];
      dcA_ijk[2][0] = gj[2];
      dcA_ijk[0][1] = gk[0];
      dcA_ijk[1][1] = gk[1];
      dcA_ijk[2][1] = gk[2];
      gfactor1 = tl_ijk.G;
      gprime1 = tl_ijk.dG;
    } else {
      angle(r_ji, dis_ji, r_jk, dis_jk, cosAng_ijk,
            dCosAngj, dCosAngk);
      int param_ijk = elem3param[itype][jtype][ktype];
      tripletParameters[param_ijk].value(cosAng_ijk, gfactor1, 1, gprime1, 1);
      dcA_ijk[0][0] = dCosAngj[0];
      dcA_ijk[1][0] = dCosAngj[1];
      dcA_ijk[2][0] = dCosAngj[2];
      dcA_ijk[0][1] = dCosAngk[0];
      dcA_ijk[1][1] = dCosAngk[1];
      dcA_ijk[2][1] = dCosAngk[2];
    }
    gfactorsq = gfactor1*gfactor1;
    gsqprime = 2.0*gfactor1*gprime1;

    // BB is Eq. 34 (a) or Eq. 10 (c) for the j atom
    // 1st DD is Eq. 11 (c) for j atom where i & k=neighbor of j

    BB += gfactorsq * betaS_jk * betaS_jk;

    // agpdpr1 is derivative of BB w.r.t. r_jk
    // app1 is derivative of BB w.r.t. cos(theta_ijk)

    agpdpr1 = 2.0 * gfactorsq * betaS_jk * dBetaS_jk / r_jk;
    app1 = betaS_jk * betaS_jk * gsqprime;
    bt_sg[nb_ij].dBB[0] -= app1*dcA_ijk[0][0];
    bt_sg[nb_ij].dBB[1] -= app1*dcA_ijk[1][0];
    bt_sg[nb_ij].dBB[2] -= app1*dcA_ijk[2][0];
    bt_sg[nb_jk].dBB[0] += app1*dcA_ijk[0][1] + agpdpr1*dis_jk[0];
    bt_sg[nb_jk].dBB[1] += app1*dcA_ijk[1][1] + agpdpr1*dis_jk[1];
    bt_sg[nb_jk].dBB[2] += app1*dcA_ijk[2][1] + agpdpr1*dis_jk[2];
  }

  AAC = AA + BB;
  for (loop = 0; loop < nb_t; loop++) {
    bt_sg[loop].dAAC[0] = bt_sg[loop].dAA[0] + bt_sg[loop].dBB[0];
    bt_sg[loop].dAAC[1] = bt_sg[loop].dAA[1] + bt_sg[loop].dBB[1];
    bt_sg[loop].dAAC[2] = bt_sg[loop].dAA[2] + bt_sg[loop].dBB[2];
  }
  bndtmp = FF + sigma_delta[param_ij]*sigma_delta[param_ij] +
    sigma_c[param_ij]*AAC + small4;
  bndtmp0 = 1.0/sqrt(bndtmp);
  sigB1 = betaS_ij*bndtmp0;

  // bndtmp1 is derivative of SigB1 w.r.t. AAC
  // bndtmp2 is derivative of SigB1 w.r.t. r_ij

  bndtmp = -0.5*bndtmp0*bndtmp0*bndtmp0;
  bndtmp1 = (bndtmp0+betaS_ij*bndtmp*2.0*betaS_ij) * dBetaS_ij/r_ij;
  bndtmp2 = betaS_ij*bndtmp*sigma_c[param_ij];
  for (loop = 0; loop < nb_t; loop++) {
    temp_kk = bt_sg[loop].temp;
    bt_sg[loop].dSigB1[0] = bndtmp2*bt_sg[loop].dAAC[0];
    bt_sg[loop].dSigB1[1] = bndtmp2*bt_sg[loop].dAAC[1];
    bt_sg[loop].dSigB1[2] = bndtmp2*bt_sg[loop].dAAC[2];
    if (temp_kk == temp_ij) {
      bt_sg[loop].dSigB1[0] += bndtmp1*dis_ij[0];
      bt_sg[loop].dSigB1[1] += bndtmp1*dis_ij[1];
      bt_sg[loop].dSigB1[2] += bndtmp1*dis_ij[2];
    }
  }

  // sigB is the final expression for (a) Eq. 6 and (b) Eq. 11

  (p_ij.bo)->value(sigB1, sigB, 1, dsigB, 1);
  if (sigma_f[param_ij] != 0.5 && sigma_k[param_ij] != 0.0) {
    part0 = (FF+0.5*AAC+small5);
    part1 = (sigma_f[param_ij]-0.5)*sigma_k[param_ij];
    part2 = 1.0-part1*EE1/part0;
    part3 = sigB*part1/part0;
    part4 = part3/part0*EE1;
    sigB *= part2;
  }

  pp1 = 2.0*betaS_ij;
  for (loop = 0; loop < nb_t; loop++) {
    bt_i = bt_sg[loop].i;
    bt_j = bt_sg[loop].j;
    xtmp[0] = x[bt_i][0]-x[bt_j][0];
    xtmp[1] = x[bt_i][1]-x[bt_j][1];
    xtmp[2] = x[bt_i][2]-x[bt_j][2];
    if (sigma_f[param_ij] == 0.5 || sigma_k[param_ij] == 0.0) {
      for (int n = 0; n < 3; n++) {
        bt_sg[loop].dSigB[n] = dsigB*bt_sg[loop].dSigB1[n];
      }
      for (int n = 0; n < 3; n++) {
        ftmp[n] = pp1*bt_sg[loop].dSigB[n];
        f[bt_i][n] -= ftmp[n];
        f[bt_j][n] += ftmp[n];
      }
      if (evflag) ev_tally_xyz(bt_i,bt_j,nlocal,newton_pair,0.0,0.0,
                               -ftmp[0],-ftmp[1],-ftmp[2],xtmp[0],xtmp[1],xtmp[2]);
    } else {
      for (int n = 0; n < 3; n++) {
        bt_sg[loop].dSigB[n] = dsigB*part2*bt_sg[loop].dSigB1[n] -
          part3*bt_sg[loop].dEE1[n] + part4*(bt_sg[loop].dFF[n] +
                                             0.5*bt_sg[loop].dAAC[n]);
      }
      for (int n = 0; n < 3; n++) {
        ftmp[n] = pp1*bt_sg[loop].dSigB[n];
        f[bt_i][n] -= ftmp[n];
        f[bt_j][n] += ftmp[n];
      }
      if (evflag) ev_tally_xyz(bt_i,bt_j,nlocal,newton_pair,0.0,0.0,
                               -ftmp[0],-ftmp[1],-ftmp[2],xtmp[0],xtmp[1],xtmp[2]);
    }
  }
  return(sigB);
}

/* ---------------------------------------------------------------------- */

/*  The formulation differs slightly to avoid negative square roots
    in the calculation of Theta_pi,ij of (a) Eq. 36 and (b) Eq. 18
    see (d) */

/* ---------------------------------------------------------------------- */

double PairBOP::PiBo(int itmp, int jtmp)
{
  double piB, ftmp[3], xtmp[3];
  int i, j, k, kp;
  int itype, jtype;
  int nb_t, nb_ij, nb_ik, nb_ikp, nb_jk, nb_jkp, bt_i, bt_j;
  int n_ji, n_jik, n_jikp, n_kikp, n_ijk, n_ijkp, n_kjkp;
  int temp_ij, temp_ik, temp_ikp, temp_jk, temp_jkp, temp_kk, temp_jik,
    temp_jikp, temp_kikp, temp_ijk, temp_ijkp, temp_kjkp;
  int *ilist, *jlist;
  int nlisti, nlistj;
  double r_ij, dis_ij[3], r_ik, dis_ik[3], r_ikp, dis_ikp[3], r_ji, dis_ji[3],
    r_jk, dis_jk[3], r_jkp, dis_jkp[3];
  double dCosAngj[3], dCosAngk[3];
  double *gj, *gk;
  double betaP_ij, dBetaP_ij, betaS_ik, dBetaS_ik, betaP_ik, dBetaP_ik,
    betaS_ikp, dBetaS_ikp, betaP_ikp, dBetaP_ikp, betaS_jk, dBetaS_jk,
    betaP_jk, dBetaP_jk, betaS_jkp, dBetaS_jkp, betaP_jkp, dBetaP_jkp;
  double cosAng_jik, dcA_jik[3][2], cosAng_jikp, dcA_jikp[3][2],
    cosAng_kikp, dcA_kikp[3][2], cosAng_ijk, dcA_ijk[3][2], cosAng_ijkp,
    dcA_ijkp[3][2], cosAng_kjkp, dcA_kjkp[3][2];
  double AA, BB, AB1, AB2, CC, BBrt, BBrtR, ABrtR1, ABrtR2, dPiB1, dPiB2,
    dPiB3, pp2;
  double cosSq, cosSq1, sinFactor, cosFactor, betaCapSq1,
    dbetaCapSq1, betaCapSq2, dbetaCapSq2, agpdpr1, agpdpr2, agpdpr3,
    app1, app2, app3, angFactor, angFactor1, angFactor2, angFactor3,
    angFactor4, angRfactor, betaCapSum, dotV, dAngR1, dAngR2;
  int loop, temp_loop, nei_loop, nei;
  int ktmp, ltmp;

  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int *iilist = list->ilist;
  int **firstneigh = list->firstneigh;

  // Loop over all local atoms for i

  piB = 0;
  if (itmp < nlocal) {
    i = iilist[itmp];
  } else {
    i = itmp;
  }

  nb_t = 0;
  memory_pi(nb_t);
  initial_pi(nb_t);

  itype = map[type[i]];
  ilist = firstneigh[i];
  nlisti = BOP_total[i];
  temp_ij = BOP_index[i]+jtmp;
  PairList1 & pl_ij = pairlist1[temp_ij];
  j = ilist[neigh_index[temp_ij]];
  jtype = map[type[j]];
  jlist = firstneigh[j];
  nlistj = BOP_total[j];
  int param_ij = elem2param[itype][jtype];
  nb_ij = nb_t;
  bt_pi[nb_ij].i = i;
  bt_pi[nb_ij].j = j;
  bt_pi[nb_ij].temp = temp_ij;
  nb_t++;
  memory_pi(nb_t);
  initial_pi(nb_t);

  for (loop = 0; loop < nlistj; loop++) {
    temp_loop = BOP_index[j] + loop;
    nei_loop = neigh_index[temp_loop];
    nei = jlist[nei_loop];
    if (x[nei][0]==x[i][0] && x[nei][1]==x[i][1] && x[nei][2]==x[i][2]) {
      n_ji = loop;
      break;
    }
  }

  dis_ij[0] = pl_ij.dis[0];
  dis_ij[1] = pl_ij.dis[1];
  dis_ij[2] = pl_ij.dis[2];
  r_ij = pl_ij.r;
  betaP_ij = pl_ij.betaP;
  dBetaP_ij = pl_ij.dBetaP;
  dis_ji[0] = -dis_ij[0];
  dis_ji[1] = -dis_ij[1];
  dis_ji[2] = -dis_ij[2];
  r_ji = r_ij;

  AA = 0.0;
  BB = 0.0;

  // if (betaP_ij * betaP_ij <= 0.000001) return(piB);

  for (ktmp = 0; ktmp < nlisti; ktmp++) {
    if (ktmp == jtmp) continue;
    temp_ik = BOP_index[i] + ktmp;
    PairList1 & pl_ik = pairlist1[temp_ik];
    k = ilist[neigh_index[temp_ik]];

    dis_ik[0] = pl_ik.dis[0];
    dis_ik[1] = pl_ik.dis[1];
    dis_ik[2] = pl_ik.dis[2];
    r_ik = pl_ik.r;
    betaS_ik = pl_ik.betaS;
    dBetaS_ik = pl_ik.dBetaS;
    betaP_ik = pl_ik.betaP;
    dBetaP_ik = pl_ik.dBetaP;

    nb_ik = nb_t;
    bt_pi[nb_ik].i = i;
    bt_pi[nb_ik].j = k;
    bt_pi[nb_ik].temp = temp_ik;
    nb_t++;
    memory_pi(nb_t);
    initial_pi(nb_t);

    if (!otfly) {
      if (jtmp < ktmp) {
        n_jik = jtmp*(2*nlisti-jtmp-1)/2 + (ktmp-jtmp)-1;
      } else {
        n_jik = ktmp*(2*nlisti-ktmp-1)/2 + (jtmp-ktmp)-1;
      }
      temp_jik = cos_index[i] + n_jik;
      TripleList & tl_jik = triplelist[temp_jik];
      if (jtmp < ktmp) {
        gj = tl_jik.dCosAngj;
        gk = tl_jik.dCosAngk;
      } else {
        gj = tl_jik.dCosAngk;
        gk = tl_jik.dCosAngj;
      }
      cosAng_jik = tl_jik.cosAng;
      dcA_jik[0][0] = gj[0];
      dcA_jik[1][0] = gj[1];
      dcA_jik[2][0] = gj[2];
      dcA_jik[0][1] = gk[0];
      dcA_jik[1][1] = gk[1];
      dcA_jik[2][1] = gk[2];
    } else {
      angle(r_ij, dis_ij, r_ik, dis_ik, cosAng_jik,
            dCosAngj, dCosAngk);
      dcA_jik[0][0] = dCosAngj[0];
      dcA_jik[1][0] = dCosAngj[1];
      dcA_jik[2][0] = dCosAngj[2];
      dcA_jik[0][1] = dCosAngk[0];
      dcA_jik[1][1] = dCosAngk[1];
      dcA_jik[2][1] = dCosAngk[2];
    }

    cosSq = cosAng_jik*cosAng_jik;
    sinFactor = 0.5*(1.0-cosSq)*pi_p[itype]*betaS_ik;
    cosFactor = 0.5*(1.0+cosSq)*betaP_ik;
    betaCapSq1 = pi_p[itype]*betaS_ik*betaS_ik-betaP_ik*betaP_ik;
    dbetaCapSq1 = 2.0*pi_p[itype]*betaS_ik*dBetaS_ik -
      2.0*betaP_ik*dBetaP_ik;

    // AA is Eq. 37 (a) and Eq. 19 (b) or i atoms
    // 1st BB is first term of Eq. 38 (a) where j and k =neighbors i

    AA = AA + sinFactor*betaS_ik + cosFactor*betaP_ik;
    BB = BB + 0.25*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*betaCapSq1;

    // agpdpr1 is derivative of AA w.r.t. for atom i w.r.t. r_ik
    // agpdpr2 is derivative of BB w.r.t. for atom i w.r.t. r_ik
    // app1 is derivative of AA w.r.t. for atom i w.r.t. cos(theta_jik)
    // app2 is derivative of BB w.r.t. for atom i w.r.t. cos(theta_jik)

    agpdpr1 = (2.0*sinFactor*dBetaS_ik+2.0*cosFactor*dBetaP_ik) / r_ik;
    agpdpr2 = 0.5*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*dbetaCapSq1 / r_ik;
    app1 = cosAng_jik*(betaP_ik*betaP_ik-pi_p[itype]*betaS_ik*betaS_ik);
    app2 = (cosSq-1.0)*cosAng_jik*betaCapSq1*betaCapSq1;
    bt_pi[nb_ij].dAA[0] += app1*dcA_jik[0][0];
    bt_pi[nb_ij].dAA[1] += app1*dcA_jik[1][0];
    bt_pi[nb_ij].dAA[2] += app1*dcA_jik[2][0];
    bt_pi[nb_ij].dBB[0] += app2*dcA_jik[0][0];
    bt_pi[nb_ij].dBB[1] += app2*dcA_jik[1][0];
    bt_pi[nb_ij].dBB[2] += app2*dcA_jik[2][0];
    bt_pi[nb_ik].dAA[0] += agpdpr1*dis_ik[0] + app1*dcA_jik[0][1];
    bt_pi[nb_ik].dAA[1] += agpdpr1*dis_ik[1] + app1*dcA_jik[1][1];
    bt_pi[nb_ik].dAA[2] += agpdpr1*dis_ik[2] + app1*dcA_jik[2][1];
    bt_pi[nb_ik].dBB[0] += agpdpr2*dis_ik[0] + app2*dcA_jik[0][1];
    bt_pi[nb_ik].dBB[1] += agpdpr2*dis_ik[1] + app2*dcA_jik[1][1];
    bt_pi[nb_ik].dBB[2] += agpdpr2*dis_ik[2] + app2*dcA_jik[2][1];

    // j and k and k' are different neighbors of i

    for (ltmp = 0; ltmp < ktmp; ltmp++) {
      if (ltmp == jtmp) continue;
      temp_ikp = BOP_index[i] + ltmp;
      PairList1 & pl_ikp = pairlist1[temp_ikp];
      kp = ilist[neigh_index[temp_ikp]];
      dis_ikp[0] = pl_ikp.dis[0];
      dis_ikp[1] = pl_ikp.dis[1];
      dis_ikp[2] = pl_ikp.dis[2];
      r_ikp =  pl_ikp.r;
      betaS_ikp = pl_ikp.betaS;
      dBetaS_ikp = pl_ikp.dBetaS;
      betaP_ikp = pl_ikp.betaP;
      dBetaP_ikp = pl_ikp.dBetaP;

      nb_ikp = nb_t;
      bt_pi[nb_ikp].i = i;
      bt_pi[nb_ikp].j = kp;
      bt_pi[nb_ikp].temp = temp_ikp;
      nb_t++;
      memory_pi(nb_t);
      initial_pi(nb_t);

      if (!otfly) {
        n_kikp = ltmp*(2*nlisti-ltmp-1)/2 + (ktmp-ltmp)-1;
        temp_kikp = cos_index[i] + n_kikp;
        TripleList & tl_kikp = triplelist[temp_kikp];
        gj = tl_kikp.dCosAngk;
        gk = tl_kikp.dCosAngj;
        cosAng_kikp = tl_kikp.cosAng;
        dcA_kikp[0][0] = gj[0];
        dcA_kikp[1][0] = gj[1];
        dcA_kikp[2][0] = gj[2];
        dcA_kikp[0][1] = gk[0];
        dcA_kikp[1][1] = gk[1];
        dcA_kikp[2][1] = gk[2];
        if (jtmp < ltmp) {
          n_jikp = jtmp*(2*nlisti-jtmp-1)/2 + (ltmp-jtmp)-1;
        } else {
          n_jikp = ltmp*(2*nlisti-ltmp-1)/2 + (jtmp-ltmp)-1;
        }
        temp_jikp = cos_index[i] + n_jikp;
        TripleList & tl_jikp = triplelist[temp_jikp];
        if (jtmp < ltmp) {
          gj = tl_jikp.dCosAngj;
          gk = tl_jikp.dCosAngk;
        } else {
          gj = tl_jikp.dCosAngk;
          gk = tl_jikp.dCosAngj;
        }
        cosAng_jikp = tl_jikp.cosAng;
        dcA_jikp[0][0] = gj[0];
        dcA_jikp[1][0] = gj[1];
        dcA_jikp[2][0] = gj[2];
        dcA_jikp[0][1] = gk[0];
        dcA_jikp[1][1] = gk[1];
        dcA_jikp[2][1] = gk[2];
      } else {
        angle(r_ik, dis_ik, r_ikp, dis_ikp, cosAng_kikp,
              dCosAngj, dCosAngk);
        dcA_kikp[0][0] = dCosAngj[0];
        dcA_kikp[1][0] = dCosAngj[1];
        dcA_kikp[2][0] = dCosAngj[2];
        dcA_kikp[0][1] = dCosAngk[0];
        dcA_kikp[1][1] = dCosAngk[1];
        dcA_kikp[2][1] = dCosAngk[2];
        angle(r_ij, dis_ij, r_ikp, dis_ikp, cosAng_jikp,
              dCosAngj, dCosAngk);
        dcA_jikp[0][0] = dCosAngj[0];
        dcA_jikp[1][0] = dCosAngj[1];
        dcA_jikp[2][0] = dCosAngj[2];
        dcA_jikp[0][1] = dCosAngk[0];
        dcA_jikp[1][1] = dCosAngk[1];
        dcA_jikp[2][1] = dCosAngk[2];
      }

      betaCapSq2 = pi_p[itype]*betaS_ikp*betaS_ikp - betaP_ikp*betaP_ikp;
      dbetaCapSq2 = 2.0*pi_p[itype]*betaS_ikp*dBetaS_ikp -
        2.0*betaP_ikp*dBetaP_ikp;
      cosSq1 = cosAng_jikp*cosAng_jikp;
      angFactor = cosAng_kikp-cosAng_jikp*cosAng_jik;
      angFactor1 = 4.0*angFactor;
      angFactor2 = -angFactor1*cosAng_jikp + 2.0*cosAng_jik*(1.0-cosSq1);
      angFactor3 = -angFactor1*cosAng_jik + 2.0*cosAng_jikp*(1.0-cosSq);
      angFactor4 = 2.0*angFactor*angFactor - (1.0-cosSq)*(1.0-cosSq1);
      betaCapSum = 0.5*betaCapSq1*betaCapSq2;

      // 2nd BB is third term of Eq. 38 (a) where j , k and k'=neighbors i

      BB += betaCapSum*angFactor4;

      // agpdpr1 is derivative of BB w.r.t. for atom i w.r.t. r_ik
      // agpdpr2 is derivative of BB w.r.t. for atom i w.r.t. r_ik'
      // app1 is derivative of BB 3rd term w.r.t. cos(theta_kik')
      // app2 is derivative of BB 3rd term w.r.t. cos(theta_jik)
      // app3 is derivative of BB 3rd term w.r.t. cos(theta_jik')

      agpdpr1 = 0.5*angFactor4*dbetaCapSq1*betaCapSq2/r_ik;
      agpdpr2 = 0.5*angFactor4*betaCapSq1*dbetaCapSq2/r_ikp;
      app1 = betaCapSum*angFactor1;
      app2 = betaCapSum*angFactor2;
      app3 = betaCapSum*angFactor3;
      bt_pi[nb_ij].dBB[0] += app2*dcA_jik[0][0] + app3*dcA_jikp[0][0];
      bt_pi[nb_ij].dBB[1] += app2*dcA_jik[1][0] + app3*dcA_jikp[1][0];
      bt_pi[nb_ij].dBB[2] += app2*dcA_jik[2][0] + app3*dcA_jikp[2][0];
      bt_pi[nb_ik].dBB[0] += agpdpr1*dis_ik[0] + app1*dcA_kikp[0][0] +
        app2*dcA_jik[0][1];
      bt_pi[nb_ik].dBB[1] += agpdpr1*dis_ik[1] + app1*dcA_kikp[1][0] +
        app2*dcA_jik[1][1];
      bt_pi[nb_ik].dBB[2] += agpdpr1*dis_ik[2] + app1*dcA_kikp[2][0] +
        app2*dcA_jik[2][1];
      bt_pi[nb_ikp].dBB[0] += agpdpr2*dis_ikp[0] + app1*dcA_kikp[0][1] +
        app3*dcA_jikp[0][1];
      bt_pi[nb_ikp].dBB[1] += agpdpr2*dis_ikp[1] + app1*dcA_kikp[1][1] +
        app3*dcA_jikp[1][1];
      bt_pi[nb_ikp].dBB[2] += agpdpr2*dis_ikp[2] + app1*dcA_kikp[2][1] +
        app3*dcA_jikp[2][1];
    }
  }

  for (ktmp = 0; ktmp < nlistj; ktmp++) {
    if (ktmp == n_ji) continue;
    temp_jk = BOP_index[j] + ktmp;
    PairList1 & pl_jk = pairlist1[temp_jk];
    k = jlist[neigh_index[temp_jk]];

    dis_jk[0] = pl_jk.dis[0];
    dis_jk[1] = pl_jk.dis[1];
    dis_jk[2] = pl_jk.dis[2];
    r_jk = pl_jk.r;
    betaS_jk = pl_jk.betaS;
    dBetaS_jk = pl_jk.dBetaS;
    betaP_jk = pl_jk.betaP;
    dBetaP_jk = pl_jk.dBetaP;

    nb_jk = nb_t;
    bt_pi[nb_jk].i = j;
    bt_pi[nb_jk].j = k;
    bt_pi[nb_jk].temp = temp_jk;
    nb_t++;
    memory_pi(nb_t);
    initial_pi(nb_t);

    if (!otfly) {
      if (n_ji < ktmp) {
        n_ijk = n_ji*(2*nlistj-n_ji-1)/2 + (ktmp-n_ji)-1;
      } else {
        n_ijk = ktmp*(2*nlistj-ktmp-1)/2 + (n_ji-ktmp)-1;
      }
      temp_ijk = cos_index[j] + n_ijk;
      TripleList & tl_ijk = triplelist[temp_ijk];
      if (n_ji < ktmp) {
        gj = tl_ijk.dCosAngj;
        gk = tl_ijk.dCosAngk;
      } else {
        gj = tl_ijk.dCosAngk;
        gk = tl_ijk.dCosAngj;
      }
      cosAng_ijk = tl_ijk.cosAng;
      dcA_ijk[0][0] = gj[0];
      dcA_ijk[1][0] = gj[1];
      dcA_ijk[2][0] = gj[2];
      dcA_ijk[0][1] = gk[0];
      dcA_ijk[1][1] = gk[1];
      dcA_ijk[2][1] = gk[2];
    } else {
      angle(r_ji, dis_ji, r_jk, dis_jk, cosAng_ijk,
            dCosAngj, dCosAngk);
      dcA_ijk[0][0] = dCosAngj[0];
      dcA_ijk[1][0] = dCosAngj[1];
      dcA_ijk[2][0] = dCosAngj[2];
      dcA_ijk[0][1] = dCosAngk[0];
      dcA_ijk[1][1] = dCosAngk[1];
      dcA_ijk[2][1] = dCosAngk[2];
    }

    cosSq = cosAng_ijk*cosAng_ijk;
    sinFactor = 0.5*(1.0-cosSq)*pi_p[jtype]*betaS_jk;
    cosFactor = 0.5*(1.0+cosSq)*betaP_jk;
    betaCapSq1 = pi_p[jtype]*betaS_jk*betaS_jk - betaP_jk*betaP_jk;
    dbetaCapSq1 = 2.0*pi_p[jtype]*betaS_jk*dBetaS_jk -
      2.0*betaP_jk*dBetaP_jk;

    // AA is Eq. 37 (a) and Eq. 19 (b) for j atoms
    // 3rd BB is 2nd term of Eq. 38 (a) where i and k =neighbors j

    AA += sinFactor*betaS_jk + cosFactor*betaP_jk;
    BB += 0.25*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*betaCapSq1;

    // agpdpr1 is derivative of AA for atom j w.r.t. r_jk
    // agpdpr2 is derivative of BB for atom j w.r.t. r_jk
    // app1 is derivative of AA for j atom w.r.t. cos(theta_ijk)
    // app2 is derivative of BB 2nd term w.r.t. cos(theta_ijk)

    agpdpr1 = (2.0*sinFactor*dBetaS_jk+2.0*cosFactor*dBetaP_jk) / r_jk;
    agpdpr2 = 0.5*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*dbetaCapSq1 / r_jk;
    app1 = cosAng_ijk*(betaP_jk*betaP_jk-pi_p[jtype]*betaS_jk*betaS_jk);
    app2 = (cosSq-1.0)*cosAng_ijk*betaCapSq1*betaCapSq1;
    bt_pi[nb_ij].dAA[0] -= app1*dcA_ijk[0][0];
    bt_pi[nb_ij].dAA[1] -= app1*dcA_ijk[1][0];
    bt_pi[nb_ij].dAA[2] -= app1*dcA_ijk[2][0];
    bt_pi[nb_ij].dBB[0] -= app2*dcA_ijk[0][0];
    bt_pi[nb_ij].dBB[1] -= app2*dcA_ijk[1][0];
    bt_pi[nb_ij].dBB[2] -= app2*dcA_ijk[2][0];
    bt_pi[nb_jk].dAA[0] += agpdpr1*dis_jk[0] + app1*dcA_ijk[0][1];
    bt_pi[nb_jk].dAA[1] += agpdpr1*dis_jk[1] + app1*dcA_ijk[1][1];
    bt_pi[nb_jk].dAA[2] += agpdpr1*dis_jk[2] + app1*dcA_ijk[2][1];
    bt_pi[nb_jk].dBB[0] += agpdpr2*dis_jk[0] + app2*dcA_ijk[0][1];
    bt_pi[nb_jk].dBB[1] += agpdpr2*dis_jk[1] + app2*dcA_ijk[1][1];
    bt_pi[nb_jk].dBB[2] += agpdpr2*dis_jk[2] + app2*dcA_ijk[2][1];

    // j is a neighbor of i and k and k' are different neighbors of j not equal to i

    for (ltmp = 0; ltmp < ktmp; ltmp++) {
      if (ltmp == n_ji) continue;
      temp_jkp = BOP_index[j]+ltmp;
      PairList1 & pl_jkp = pairlist1[temp_jkp];
      kp = jlist[neigh_index[temp_jkp]];
      dis_jkp[0] = pl_jkp.dis[0];
      dis_jkp[1] = pl_jkp.dis[1];
      dis_jkp[2] = pl_jkp.dis[2];
      r_jkp = pl_jkp.r;
      betaS_jkp = pl_jkp.betaS;
      dBetaS_jkp = pl_jkp.dBetaS;
      betaP_jkp = pl_jkp.betaP;
      dBetaP_jkp = pl_jkp.dBetaP;

      nb_jkp = nb_t;
      bt_pi[nb_jkp].i = j;
      bt_pi[nb_jkp].j = kp;
      bt_pi[nb_jkp].temp = temp_jkp;
      nb_t++;
      memory_pi(nb_t);
      initial_pi(nb_t);

      if (!otfly) {
        n_kjkp = ltmp*(2*nlistj-ltmp-1)/2 + (ktmp-ltmp)-1;
        temp_kjkp = cos_index[j] + n_kjkp;
        TripleList & tl_kjkp = triplelist[temp_kjkp];
        gj = tl_kjkp.dCosAngk;
        gk = tl_kjkp.dCosAngj;
        cosAng_kjkp = tl_kjkp.cosAng;
        dcA_kjkp[0][0] = gj[0];
        dcA_kjkp[1][0] = gj[1];
        dcA_kjkp[2][0] = gj[2];
        dcA_kjkp[0][1] = gk[0];
        dcA_kjkp[1][1] = gk[1];
        dcA_kjkp[2][1] = gk[2];
        if (n_ji < ltmp) {
          n_ijkp = n_ji*(2*nlistj-n_ji-1)/2 + (ltmp-n_ji)-1;
        } else {
          n_ijkp = ltmp*(2*nlistj-ltmp-1)/2 + (n_ji-ltmp)-1;
        }
        temp_ijkp = cos_index[j] + n_ijkp;
        TripleList & tl_ijkp = triplelist[temp_ijkp];
        if (n_ji < ltmp) {
          gj = tl_ijkp.dCosAngj;
          gk = tl_ijkp.dCosAngk;
        } else {
          gj = tl_ijkp.dCosAngk;
          gk = tl_ijkp.dCosAngj;
        }
        cosAng_ijkp = tl_ijkp.cosAng;
        dcA_ijkp[0][0] = gj[0];
        dcA_ijkp[1][0] = gj[1];
        dcA_ijkp[2][0] = gj[2];
        dcA_ijkp[0][1] = gk[0];
        dcA_ijkp[1][1] = gk[1];
        dcA_ijkp[2][1] = gk[2];
      } else {
        angle(r_jk, dis_jk, r_jkp, dis_jkp, cosAng_kjkp,
              dCosAngj, dCosAngk);
        dcA_kjkp[0][0] = dCosAngj[0];
        dcA_kjkp[1][0] = dCosAngj[1];
        dcA_kjkp[2][0] = dCosAngj[2];
        dcA_kjkp[0][1] = dCosAngk[0];
        dcA_kjkp[1][1] = dCosAngk[1];
        dcA_kjkp[2][1] = dCosAngk[2];
        angle(r_ji, dis_ji, r_jkp, dis_jkp, cosAng_ijkp,
              dCosAngj, dCosAngk);
        dcA_ijkp[0][0] = dCosAngj[0];
        dcA_ijkp[1][0] = dCosAngj[1];
        dcA_ijkp[2][0] = dCosAngj[2];
        dcA_ijkp[0][1] = dCosAngk[0];
        dcA_ijkp[1][1] = dCosAngk[1];
        dcA_ijkp[2][1] = dCosAngk[2];
      }

      betaCapSq2 = pi_p[jtype]*betaS_jkp*betaS_jkp -betaP_jkp*betaP_jkp;
      dbetaCapSq2 = 2.0*pi_p[jtype]*betaS_jkp*dBetaS_jkp -
        2.0*betaP_jkp*dBetaP_jkp;
      cosSq1 = cosAng_ijkp*cosAng_ijkp;
      angFactor = cosAng_kjkp-cosAng_ijkp*cosAng_ijk;
      angFactor1 = 4.0*angFactor;
      angFactor2 = -angFactor1*cosAng_ijkp + 2.0*cosAng_ijk*(1.0-cosSq1);
      angFactor3 = -angFactor1*cosAng_ijk + 2.0*cosAng_ijkp*(1.0-cosSq);
      angFactor4 = 2.0*angFactor*angFactor-(1.0-cosSq)*(1.0-cosSq1);
      betaCapSum = 0.5*betaCapSq1*betaCapSq2;

      // 4th BB is 4th term of Eq. 38 (a) where i , k and k' =neighbors j

      BB += betaCapSum*angFactor4;

      // app1 is derivative of BB 4th term w.r.t. cos(theta_kjk')
      // app2 is derivative of BB 4th term w.r.t. cos(theta_ijk)
      // app3 is derivative of BB 4th term w.r.t. cos(theta_ijk')
      // agpdpr1 is derivative of BB 4th term for atom j w.r.t. r_jk
      // agpdpr2 is derivative of BB 4th term for atom j w.r.t. r_jk'

      agpdpr1 = 0.5*angFactor4*dbetaCapSq1*betaCapSq2 / r_jk;
      agpdpr2 = 0.5*angFactor4*betaCapSq1*dbetaCapSq2 / r_jkp;
      app1 = betaCapSum*angFactor1;
      app2 = betaCapSum*angFactor2;
      app3 = betaCapSum*angFactor3;
      bt_pi[nb_ij].dBB[0] -= app3*dcA_ijkp[0][0] + app2*dcA_ijk[0][0];
      bt_pi[nb_ij].dBB[1] -= app3*dcA_ijkp[1][0] + app2*dcA_ijk[1][0];
      bt_pi[nb_ij].dBB[2] -= app3*dcA_ijkp[2][0] + app2*dcA_ijk[2][0];
      bt_pi[nb_jk].dBB[0] += agpdpr1*dis_jk[0] + app1*dcA_kjkp[0][0] +
        app2*dcA_ijk[0][1];
      bt_pi[nb_jk].dBB[1] += agpdpr1*dis_jk[1] + app1*dcA_kjkp[1][0] +
        app2*dcA_ijk[1][1];
      bt_pi[nb_jk].dBB[2] += agpdpr1*dis_jk[2] + app1*dcA_kjkp[2][0] +
        app2*dcA_ijk[2][1];
      bt_pi[nb_jkp].dBB[0] += agpdpr2*dis_jkp[0] + app1*dcA_kjkp[0][1] +
        app3*dcA_ijkp[0][1];
      bt_pi[nb_jkp].dBB[1] += agpdpr2*dis_jkp[1] + app1*dcA_kjkp[1][1] +
        app3*dcA_ijkp[1][1];
      bt_pi[nb_jkp].dBB[2] += agpdpr2*dis_jkp[2] + app1*dcA_kjkp[2][1] +
        app3*dcA_ijkp[2][1];
    }

    // j and k' are different neighbors of i and k is a neighbor of j not equal to i

    for (ltmp = 0; ltmp < nlisti; ltmp++) {
      if (ltmp == jtmp) continue;
      temp_ikp = BOP_index[i] + ltmp;
      PairList1 & pl_ikp = pairlist1[temp_ikp];
      kp=ilist[neigh_index[temp_ikp]];
      dis_ikp[0] = pl_ikp.dis[0];
      dis_ikp[1] = pl_ikp.dis[1];
      dis_ikp[2] = pl_ikp.dis[2];
      r_ikp = pl_ikp.r;
      betaS_ikp = pl_ikp.betaS;
      dBetaS_ikp = pl_ikp.dBetaS;
      betaP_ikp = pl_ikp.betaP;
      dBetaP_ikp = pl_ikp.dBetaP;

      nb_ikp = nb_t;
      bt_pi[nb_ikp].i = i;
      bt_pi[nb_ikp].j = kp;
      bt_pi[nb_ikp].temp = temp_ikp;
      nb_t++;
      memory_pi(nb_t);
      initial_pi(nb_t);

      if (!otfly) {
        if (jtmp < ltmp) {
          n_jikp = jtmp*(2*nlisti-jtmp-1)/2 + (ltmp-jtmp)-1;
        } else {
          n_jikp = ltmp*(2*nlisti-ltmp-1)/2 + (jtmp-ltmp)-1;
        }
        temp_jikp = cos_index[i] + n_jikp;
        TripleList & tl_jikp = triplelist[temp_jikp];
        if (jtmp < ltmp) {
          gj = tl_jikp.dCosAngj;
          gk = tl_jikp.dCosAngk;
        } else {
          gj = tl_jikp.dCosAngk;
          gk = tl_jikp.dCosAngj;
        }
        cosAng_jikp = tl_jikp.cosAng;
        dcA_jikp[0][0] = gj[0];
        dcA_jikp[1][0] = gj[1];
        dcA_jikp[2][0] = gj[2];
        dcA_jikp[0][1] = gk[0];
        dcA_jikp[1][1] = gk[1];
        dcA_jikp[2][1] = gk[2];
      } else {
        angle(r_ij, dis_ij, r_ikp, dis_ikp, cosAng_jikp,
              dCosAngj, dCosAngk);
        dcA_jikp[0][0] = dCosAngj[0];
        dcA_jikp[1][0] = dCosAngj[1];
        dcA_jikp[2][0] = dCosAngj[2];
        dcA_jikp[0][1] = dCosAngk[0];
        dcA_jikp[1][1] = dCosAngk[1];
        dcA_jikp[2][1] = dCosAngk[2];
      }

      betaCapSq2 = pi_p[itype]*betaS_ikp*betaS_ikp - betaP_ikp*betaP_ikp;
      dbetaCapSq2 = 2.0*pi_p[itype]*betaS_ikp*dBetaS_ikp -
        2.0*betaP_ikp*dBetaP_ikp;
      dotV = (dis_jk[0]*dis_ikp[0] + dis_jk[1]*dis_ikp[1] +
              dis_jk[2]*dis_ikp[2]) / (r_jk*r_ikp);
      cosSq1 = cosAng_jikp*cosAng_jikp;
      angFactor = dotV + cosAng_jikp*cosAng_ijk;
      angRfactor = 4.0*angFactor*dotV;
      dAngR1 = -angRfactor/r_jk;
      dAngR2 = -angRfactor/r_ikp;
      angFactor1 = 4.0*angFactor*cosAng_jikp + 2.0*cosAng_ijk*(1.0-cosSq1);
      angFactor2 = 4.0*angFactor*cosAng_ijk + 2.0*cosAng_jikp*(1.0-cosSq);
      angFactor3 = 2.0*angFactor*angFactor - (1.0-cosSq)*(1.0-cosSq1);
      betaCapSum = 0.5*betaCapSq1*betaCapSq2;

      // 5th BB is 5th term of Eq. 38 (a) Eq. 21 (b) where i , k and k' =neighbors j

      BB += betaCapSum*angFactor3;

      // app1 is derivative of BB 5th term w.r.t. cos(theta_jik)
      // app2 is derivative of BB 5th term w.r.t. cos(theta_jik')
      // agpdpr1 is derivative of BB 5th term for atom j w.r.t. r_jk
      // agpdpr2 is derivative of BB 5th term for atom j w.r.t. r_ik'
      // agpdpr3 is derivative of BB 5th term for atom j w.r.t. dot(r_ik',r_jk)

      agpdpr1 = (0.5*angFactor3*dbetaCapSq1*betaCapSq2 +
                 betaCapSum*dAngR1) / r_jk;
      agpdpr2 = (0.5*angFactor3*betaCapSq1*dbetaCapSq2 +
                 betaCapSum*dAngR2) / r_ikp;
      agpdpr3 = 4.0*betaCapSum*angFactor/(r_ikp*r_jk);
      app1 = betaCapSum*angFactor1;
      app2 = betaCapSum*angFactor2;
      bt_pi[nb_ij].dBB[0] += app2*dcA_jikp[0][0] - app1*dcA_ijk[0][0];
      bt_pi[nb_ij].dBB[1] += app2*dcA_jikp[1][0] - app1*dcA_ijk[1][0];
      bt_pi[nb_ij].dBB[2] += app2*dcA_jikp[2][0] - app1*dcA_ijk[2][0];
      bt_pi[nb_ikp].dBB[0]+= agpdpr2*dis_ikp[0] + agpdpr3*dis_jk[0] +
        app2*dcA_jikp[0][1];
      bt_pi[nb_ikp].dBB[1] += agpdpr2*dis_ikp[1] + agpdpr3*dis_jk[1] +
        app2*dcA_jikp[1][1];
      bt_pi[nb_ikp].dBB[2] += agpdpr2*dis_ikp[2] + agpdpr3*dis_jk[2] +
        app2*dcA_jikp[2][1];
      bt_pi[nb_jk].dBB[0] += agpdpr1*dis_jk[0] + agpdpr3*dis_ikp[0] +
        app1*dcA_ijk[0][1];
      bt_pi[nb_jk].dBB[1] += agpdpr1*dis_jk[1] + agpdpr3*dis_ikp[1] +
        app1*dcA_ijk[1][1];
      bt_pi[nb_jk].dBB[2] += agpdpr1*dis_jk[2] + agpdpr3*dis_ikp[2] +
        app1*dcA_ijk[2][1];
    }
  }
  pp2 = 2.0*betaP_ij;
  CC = betaP_ij*betaP_ij + pi_delta[param_ij]*pi_delta[param_ij];
  BBrt = sqrt(BB+small6);
  AB1 = CC + pi_c[param_ij]*(AA+BBrt) + small7;
  AB2 = CC + pi_c[param_ij]*(AA-BBrt+sqrt(small6)) + small7;
  BBrtR = 1.0/BBrt;
  ABrtR1 = 1.0/sqrt(AB1);
  ABrtR2 = 1.0/sqrt(AB2);
  piB = (ABrtR1+ABrtR2)*pi_a[param_ij]*betaP_ij;

  // dPiB1 is derivative of piB w.r.t. AA
  // dPiB2 is derivative of piB w.r.t. BB
  // dPiB3 is derivative of piB w.r.t. r_ij

  dPiB1 = -0.5*(pow(ABrtR1,3)+pow(ABrtR2,3))*pi_c[param_ij]*pi_a[param_ij]*
    betaP_ij;
  dPiB2 = 0.25*BBrtR*(pow(ABrtR2,3)-pow(ABrtR1,3))*pi_c[param_ij]*
    pi_a[param_ij]*betaP_ij;
  dPiB3 = ((ABrtR1+ABrtR2)*pi_a[param_ij]-(pow(ABrtR1,3)+pow(ABrtR2,3))*
           pi_a[param_ij]*betaP_ij*betaP_ij)*dBetaP_ij / r_ij;
  for (loop = 0; loop < nb_t; loop++) {
    temp_kk = bt_pi[loop].temp;
    bt_pi[loop].dPiB[0] = dPiB1*bt_pi[loop].dAA[0] + dPiB2*bt_pi[loop].dBB[0];
    bt_pi[loop].dPiB[1] = dPiB1*bt_pi[loop].dAA[1] + dPiB2*bt_pi[loop].dBB[1];
    bt_pi[loop].dPiB[2] = dPiB1*bt_pi[loop].dAA[2] + dPiB2*bt_pi[loop].dBB[2];
    if (temp_kk == temp_ij) {
      bt_pi[loop].dPiB[0] += dPiB3*dis_ij[0];
      bt_pi[loop].dPiB[1] += dPiB3*dis_ij[1];
      bt_pi[loop].dPiB[2] += dPiB3*dis_ij[2];
    }
  }

  for (loop = 0; loop < nb_t; loop++) {
    bt_i = bt_pi[loop].i;
    bt_j = bt_pi[loop].j;
    xtmp[0] = x[bt_i][0]-x[bt_j][0];
    xtmp[1] = x[bt_i][1]-x[bt_j][1];
    xtmp[2] = x[bt_i][2]-x[bt_j][2];
    for (int n = 0; n <3; n++) {
      ftmp[n] = pp2*bt_pi[loop].dPiB[n];
      f[bt_i][n] -= ftmp[n];
      f[bt_j][n] += ftmp[n];
    }
    if (evflag) ev_tally_xyz(bt_i,bt_j,nlocal,newton_pair,0.0,0.0,
                             -ftmp[0],-ftmp[1],-ftmp[2],xtmp[0],xtmp[1],xtmp[2]);
  }
  return(piB);
}

/* ----------------------------------------------------------------------
   read BOP potential file
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void PairBOP::read_table(char *filename)
{
  int nbuf;
  double *singletable;
  int nr, nBOt, ntheta, npower, format;
  double *rcut = nullptr;
  double ****gpara = nullptr;
  PotentialFileReader *reader = nullptr;

  if (bop_elements) {
    for (int i = 0; i < bop_types; i++) delete[] bop_elements[i];
    delete[] bop_elements;
  }
  delete[] bop_masses;

  if (comm->me == 0) {
    try {
      reader = new PotentialFileReader(lmp, filename, "BOP");
      bop_types = reader->next_int();
      if (bop_types <= 0)
        error->one(FLERR,fmt::format("BOP potential file with {} "
                                     "elements",bop_types));

      bop_elements = new char*[bop_types];
      bop_masses = new double[bop_types];
      for (int i=0; i < bop_types; ++i) {
        ValueTokenizer values = reader->next_values(3);
        values.next_int();                    //  element number (ignored)
        bop_masses[i] = values.next_double(); //  element mass
        bop_elements[i] = utils::strdup(values.next_string());
      }
    } catch (TokenizerException &e) {
      error->one(FLERR,"Error reading BOP potential file: {}",e.what());
    }
  }
  MPI_Bcast(&bop_types,1,MPI_INT,0,world);
  npairs = bop_types * (bop_types + 1) / 2;
  ntriples = bop_types * bop_types * bop_types;
  allocate();
  memory->create(rcut,npairs,"BOP:rcut");

  // copy element labels and masses to all MPI ranks for use with
  // write_tables() and to set the per-type masses
  if (comm->me != 0) {
    bop_elements = new char*[bop_types];
    bop_masses = new double[bop_types];
  }
  for (int i = 0; i < bop_types; ++i) {
    int n=0;
    if (comm->me == 0) n = strlen(bop_elements[i])+1;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (comm->me != 0) bop_elements[i] = new char[n];
    MPI_Bcast(bop_elements[i],n,MPI_CHAR,0,world);
  }
  MPI_Bcast(bop_masses, bop_types, MPI_DOUBLE, 0, world);

  if (comm->me == 0) {
    try {
      //  continue reading from already opened potential file

      //  2 or 3 values in next line
      //  3 values means either table with npower = 2 if second value > 10
      //    or power function with power set to 2 < value < = 10.
      //  2 values always means parabolic function.
      ValueTokenizer values = reader->next_values(2);
      format = values.count();

      switch (format) {
      case 3:
        nr = values.next_int();
        ntheta = values.next_int();
        nBOt = values.next_int();
        if (ntheta > 10) {
          npower = 2;
        } else {
          npower = ntheta;
          ntheta = 2000;
        }
        break;

      case 2:
        nr = values.next_int();
        ntheta = 2000;
        nBOt = values.next_int();
        npower = 2;
        break;

      default:
        error->one(FLERR,"Unsupported BOP potential file format");
      }

      //  read seven "small" values in single line
      values = reader->next_values(7);
      small1  = values.next_double();
      small2  = values.next_double();
      small3g = values.next_double();
      small4  = values.next_double();
      small5  = values.next_double();
      small6  = values.next_double();
      small7  = values.next_double();

      for (int i = 0; i < bop_types; ++i)
        pi_p[i] = reader->next_double();

      cutmax = 0.0;
      for (int i = 0; i < npairs; ++i) {
        rcut[i] = reader->next_double();
        cutmax = MAX(rcut[i],cutmax);

        values = reader->next_values(4);
        sigma_c[i] = values.next_double();
        sigma_a[i] = values.next_double();
        pi_c[i]    = values.next_double();
        pi_a[i]    = values.next_double();

        values = reader->next_values(2);
        sigma_delta[i] = values.next_double();
        pi_delta[i]    = values.next_double();

        values = reader->next_values(3);
        sigma_f[i] = values.next_double();
        sigma_k[i] = values.next_double();
        small3[i]  = values.next_double();
      }
    } catch (TokenizerException &e) {
      error->one(FLERR,"Error reading BOP potential file: {}",e.what());
    }
  }

  MPI_Bcast(&nr,1,MPI_INT,0,world);
  MPI_Bcast(&ntheta,1,MPI_INT,0,world);
  MPI_Bcast(&nBOt,1,MPI_INT,0,world);
  MPI_Bcast(&npower,1,MPI_INT,0,world);
  MPI_Bcast(&small1,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small2,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small3g,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small4,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small5,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small6,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small7,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cutmax,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&pi_p[0],bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcut[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_c[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_a[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_c[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_a[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_delta[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_delta[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_f[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_k[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&small3[0],npairs,MPI_DOUBLE,0,world);

  memory->create(gpara,bop_types,bop_types,bop_types,npower+1,"BOP:gpara");
  singletable = new double[ntheta];

  nbuf = 0;
  for (int i = 0; i < bop_types; i++) {
    for (int j = 0; j < bop_types; j++) {
      for (int k = j; k < bop_types; k++) {
        if (comm->me == 0) {
          if (format == 3 && npower <= 2) {
            reader->next_dvector(singletable, ntheta);
          } else {
            reader->next_dvector(gpara[j][i][k], npower+1);
          }
          for (int n = 0; n < ntheta; n++) {
            double arg = -1.0 + 2.0 * n / (ntheta - 1.0);
            singletable[n] = gpara[j][i][k][npower];
            for (int m = npower; m > 0; m--) {
              singletable[n] = arg * singletable[n] + gpara[j][i][k][m-1];
            }
          }
        }
        MPI_Bcast(singletable,ntheta,MPI_DOUBLE,0,world);
        tripletParameters[nbuf].set_values(ntheta, -1.0, 1.0, singletable);
        elem3param[j][i][k] = nbuf;
        if (k != j) elem3param[k][i][j] = nbuf;
        nbuf++;
      }
    }
  }
  delete[] singletable;

  singletable = new double[nr];
  for (int i = 0; i < npairs; i++) {
    PairParameters &p = pairParameters[i];
    if (comm->me == 0) reader->next_dvector(singletable, nr);
    MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
    p.rep = new TabularFunction();
    (p.rep)->set_values(nr, 0.0, rcut[i], singletable);
    p.cutB = rcut[i];
    p.cutBsq = rcut[i]*rcut[i];
  }

  for (int i = 0; i < npairs; i++) {
    PairParameters &p = pairParameters[i];
    if (comm->me == 0) reader->next_dvector(singletable, nr);
    MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
    p.betaS = new TabularFunction();
    (p.betaS)->set_values(nr, 0.0, rcut[i], singletable);
  }

  for (int i = 0; i < npairs; i++) {
    PairParameters &p = pairParameters[i];
    if (comm->me == 0) reader->next_dvector(singletable, nr);
    MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
    p.betaP = new TabularFunction();
    (p.betaP)->set_values(nr, 0.0, rcut[i], singletable);
  }
  delete[] singletable;

  singletable = new double[nBOt];
  for (int i = 0; i < npairs; i++) {
    PairParameters &p = pairParameters[i];
    if (comm->me == 0) reader->next_dvector(singletable, nBOt);
    MPI_Bcast(singletable,nBOt,MPI_DOUBLE,0,world);
    p.bo = new TabularFunction();
    (p.bo)->set_values(nBOt, 0.0, 1.0, singletable);
  }
  delete[] singletable;

  nbuf = 0;
  for (int i = 0; i < bop_types; i++) {
    for (int j = i; j <= i; j++) {
      elem2param[i][j] = nbuf;
      nbuf++;
    }
  }
  for (int i = 0; i < bop_types - 1; i++) {
    for (int j = i + 1; j < bop_types; j++) {
      elem2param[i][j] = nbuf;
      elem2param[j][i] = nbuf;
      nbuf++;
    }
  }

  if (comm->me == 0) {
    reader->next_dvector(pro_delta,bop_types);
    reader->next_dvector(pro,bop_types);
    reader->next_dvector(rcut,npairs);
  }
  MPI_Bcast(&pro_delta[0],bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&pro[0],bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcut[0],npairs,MPI_DOUBLE,0,world);

  double add_cut = rcut[0];
  for (int i = 0; i < npairs; i++) {
    if (add_cut < rcut[i]) add_cut = rcut[i];
    PairParameters &p = pairParameters[i];
    p.cutL = rcut[i];
    p.cutLsq = rcut[i]*rcut[i];
  }

  if (add_cut > 0.0) {
    singletable = new double[nr];
    for (int i = 0; i < npairs; i++) {
      PairParameters &p = pairParameters[i];
      if (comm->me == 0) reader->next_dvector(singletable,nr);
      MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
      nbuf = 0;
      for (int j = 0; j < nr; j++) {
        if (singletable[j] != 0) nbuf = 1;
      }
      if (nbuf == 0) {
        rcut[i] = 0.0;
        p.cutL = rcut[i];
        p.cutLsq = rcut[i]*rcut[i];
      }
      if (rcut[i] != 0.0) {
        p.cphi = new TabularFunction();
        (p.cphi)->set_values(nr, 0.0, rcut[i], singletable);
      }
    }
    delete[] singletable;
  }

  memory->destroy(rcut);
  memory->destroy(gpara);

  if (comm->me == 0) delete reader;

#if defined(LMP_BOP_WRITE_TABLES)
  // for debugging, call write_tables() to check the tabular functions
  if (comm->me == 1) {
    write_tables(51);
  }
#endif
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairBOP::memory_usage()
{
  return(bytes);
}

/* ---------------------------------------------------------------------- */

void PairBOP::memory_sg(int n)
{
  if (bt_sg) {
    if (sglimit <= n) {
      sglimit += 500;
      memory->grow(bt_sg,sglimit,"BOP:bt_sg");
      bytes += 500 * sizeof(B_SG);
    }
  } else {
    sglimit = 2500;
    memory->create(bt_sg,sglimit,"BOP:bt_sg");
    bytes += sglimit * sizeof(B_SG);
  }
}

/* ---------------------------------------------------------------------- */

void PairBOP::memory_pi(int n)
{
  if (bt_pi) {
    if (pilimit <= n) {
      pilimit += 500;
      memory->grow(bt_pi,pilimit,"BOP:bt_pi");
      bytes += 500 * sizeof(B_PI);
    }
  } else {
    pilimit = 2500;
    memory->create(bt_pi,pilimit,"BOP:bt_pi");
    bytes += 2500 * sizeof(B_PI);
  }
}

/* ---------------------------------------------------------------------- */

void PairBOP::initial_sg(int n)
{
  B_SG & at = bt_sg[n];
  memset(&at, 0, sizeof(struct B_SG));
  at.i = -1;
  at.j = -1;
  at.temp = -1;
}

/* ---------------------------------------------------------------------- */

void PairBOP::initial_pi(int n)
{
  B_PI & at = bt_pi[n];
  memset(&at, 0, sizeof(struct B_PI));
  at.i = -1;
  at.j = -1;
  at.temp = -1;
}

/* ---------------------------------------------------------------------- */
#if defined(LMP_BOP_WRITE_TABLES)
void PairBOP::write_tables(int npts)
{
  FILE* fp =  nullptr;
  double  xmin,xmax,x,uf,vf,wf,ufp,vfp,wfp;
  std::string filename;

  for (int i = 0; i < bop_types; i++) {
    for (int j = 0; j < bop_types; j++) {
      int param = elem2param[i][j];
      PairParameters & pair = pairParameters[param];

      filename = fmt::format("{}{}_Pair_SPR_{}",bop_elements[i],
                             bop_elements[j],comm->me);

      fp = fopen(filename.c_str(), "w");
      xmin = (pair.betaS)->get_xmin();
      xmax = (pair.betaS)->get_xmax();
      xmax = xmin + (xmax - xmin) * 1.1;
      xmin = 1.0;
      for (int k = 0; k < npts; k++) {
        x = xmin + (xmax-xmin) * k / (npts-1);
        (pair.betaS)->value(x, uf, 1, ufp, 1);
        (pair.betaP)->value(x, vf, 1, vfp, 1);
        (pair.rep)->value(x, wf, 1, wfp, 1);
        fprintf(fp,"%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f \n",
                x,uf,vf,wf,ufp,vfp,wfp);
      }
      fclose(fp);

      if (pair.cutL != 0) {
        filename = fmt::format("{}{}_Pair_L_{}",bop_elements[i],
                               bop_elements[j],comm->me);
        fp = fopen(filename.c_str(), "w");
        xmin = (pair.cphi)->get_xmin();
        xmax = (pair.cphi)->get_xmax();
        xmax = xmin + (xmax - xmin) * 1.1;
        xmin = 1.0;
        for (int k = 0; k < npts; k++) {
          x = xmin + (xmax-xmin) * k / (npts-1);
          (pair.cphi)->value(x, uf, 1, ufp, 1);
          fprintf(fp,"%12.4f %12.4f %12.4f \n",x,uf,ufp);
        }
        fclose(fp);
      }
      filename = fmt::format("{}{}_Pair_BO_{}", bop_elements[i],
                             bop_elements[j], comm->me);
      fp = fopen(filename.c_str(), "w");
      xmin = (pair.bo)->get_xmin();
      xmax = (pair.bo)->get_xmax();
      for (int k = 0; k < npts; k++) {
        x = xmin + (xmax-xmin) * k / (npts-1);
        (pair.bo)->value(x, uf, 1, ufp, 1);
        fprintf(fp,"%12.4f %12.4f %12.4f \n",x,uf,ufp);
      }
      fclose(fp);
    }
  }

  for (int i = 0; i < bop_types; i++) {
    for (int j = 0; j < bop_types; j++) {
      for (int k = 0; k < bop_types; k++) {
        filename = fmt::format("{}{}{}_Triple_G_{}", bop_elements[i],
                               bop_elements[j],bop_elements[k],comm->me);
        fp = fopen(filename.c_str(), "w");
        int param = elem3param[i][j][k];
        auto &G = tripletParameters[param];
        xmin = G.get_xmin();
        xmax = G.get_xmax();
        for (int n = 0; n < npts; n++) {
          x = xmin + (xmax-xmin) * n / (npts-1);
          G.value(x, uf, 1, ufp, 1);
          fprintf(fp,"%12.4f %12.4f %12.4f \n",x,uf,ufp);
        }
        fclose(fp);
      }
    }
  }
}
#endif

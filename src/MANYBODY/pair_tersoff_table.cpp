// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Luca Ferraro (CASPUR)
   email: luca.ferraro@caspur.it

   Tersoff Potential
   References: (referenced as tersoff_2 functional form in LAMMPS manual)
    1) Tersoff, Phys. Rev. B 39, 5566 (1988)
------------------------------------------------------------------------- */

#include "pair_tersoff_table.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4

#define GRIDSTART 0.1
#define GRIDDENSITY_FCUTOFF 5000
#define GRIDDENSITY_EXP 12000
#define GRIDDENSITY_GTETA 12000
#define GRIDDENSITY_BIJ 7500

// max number of interaction per atom for environment potential

#define leadingDimensionInteractionList 64

/* ---------------------------------------------------------------------- */

PairTersoffTable::PairTersoffTable(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  params = nullptr;
  allocated = 0;

  preGtetaFunction = preGtetaFunctionDerived = nullptr;
  preCutoffFunction = preCutoffFunctionDerived = nullptr;
  exponential = nullptr;
  gtetaFunction = nullptr;
  gtetaFunctionDerived = nullptr;
  cutoffFunction = nullptr;
  cutoffFunctionDerived = nullptr;
  betaZetaPower = nullptr;
  betaZetaPowerDerived = nullptr;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairTersoffTable::~PairTersoffTable()
{
  memory->destroy(params);
  memory->destroy(elem3param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
  deallocateGrids();
  PairTersoffTable::deallocatePreLoops();
}

/* ---------------------------------------------------------------------- */

void PairTersoffTable::compute(int eflag, int vflag)
{
  int i,j,k,ii,inum,jnum;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  double xtmp,ytmp,ztmp;
  double fxtmp,fytmp,fztmp;
  int *ilist,*jlist,*numneigh,**firstneigh;


  int interpolIDX;
  double directorCos_ij_x, directorCos_ij_y, directorCos_ij_z, directorCos_ik_x, directorCos_ik_y, directorCos_ik_z;
  double invR_ij, invR_ik, cosTeta;
  double repulsivePotential, attractivePotential;
  double exponentRepulsivePotential, exponentAttractivePotential,interpolTMP,interpolDeltaX,interpolY1;
  double interpolY2, cutoffFunctionIJ, attractiveExponential, repulsiveExponential, cutoffFunctionDerivedIJ,zeta;
  double gtetaFunctionIJK,gtetaFunctionDerivedIJK,cutoffFunctionIK;
  double cutoffFunctionDerivedIK,factor_force3_ij,factor_1_force3_ik;
  double factor_2_force3_ik,betaZetaPowerIJK,betaZetaPowerDerivedIJK,factor_force_tot;
  double factor_force_ij;
  double gtetaFunctionDerived_temp,gtetaFunction_temp;

  double evdwl = 0.0;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (jnum > leadingDimensionInteractionList)
      error->one(FLERR,"Too many neighbors for interaction list: {} vs {}.\n"
                 "Check your system or increase 'leadingDimensionInteractionList'",
                 jnum, leadingDimensionInteractionList);

    // Pre-calculate gteta and cutoff function
    for (int neighbor_j = 0; neighbor_j < jnum; neighbor_j++) {

      double dr_ij[3], r_ij;

      j = jlist[neighbor_j];
      j &= NEIGHMASK;

      dr_ij[0] = xtmp - x[j][0];
      dr_ij[1] = ytmp - x[j][1];
      dr_ij[2] = ztmp - x[j][2];
      r_ij = dr_ij[0]*dr_ij[0] + dr_ij[1]*dr_ij[1] + dr_ij[2]*dr_ij[2];

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];

      if (r_ij > params[ijparam].cutsq) continue;

      r_ij = sqrt(r_ij);

      invR_ij = 1.0 / r_ij;

      directorCos_ij_x = invR_ij * dr_ij[0];
      directorCos_ij_y = invR_ij * dr_ij[1];
      directorCos_ij_z = invR_ij * dr_ij[2];

      // preCutoffFunction
      interpolDeltaX =  r_ij - GRIDSTART;
      interpolTMP = (interpolDeltaX * GRIDDENSITY_FCUTOFF);
      interpolIDX = (int) interpolTMP;
      interpolY1 = cutoffFunction[itype][jtype][interpolIDX];
      interpolY2 = cutoffFunction[itype][jtype][interpolIDX+1];
      preCutoffFunction[neighbor_j] = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);
      // preCutoffFunctionDerived
      interpolY1 = cutoffFunctionDerived[itype][jtype][interpolIDX];
      interpolY2 = cutoffFunctionDerived[itype][jtype][interpolIDX+1];
      preCutoffFunctionDerived[neighbor_j] = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);


      for (int neighbor_k = neighbor_j + 1; neighbor_k < jnum; neighbor_k++) {
        double dr_ik[3], r_ik;

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = xtmp -x[k][0];
        dr_ik[1] = ytmp -x[k][1];
        dr_ik[2] = ztmp -x[k][2];
        r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        r_ik = sqrt(r_ik);

        invR_ik = 1.0 / r_ik;

        directorCos_ik_x = invR_ik * dr_ik[0];
        directorCos_ik_y = invR_ik * dr_ik[1];
        directorCos_ik_z = invR_ik * dr_ik[2];

        cosTeta = directorCos_ij_x * directorCos_ik_x + directorCos_ij_y * directorCos_ik_y + directorCos_ij_z * directorCos_ik_z;

        // preGtetaFunction
        interpolDeltaX=cosTeta+1.0;
        interpolTMP = (interpolDeltaX * GRIDDENSITY_GTETA);
        interpolIDX = (int) interpolTMP;
        interpolY1 = gtetaFunction[itype][interpolIDX];
        interpolY2 = gtetaFunction[itype][interpolIDX+1];
        gtetaFunction_temp = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);
        // preGtetaFunctionDerived
        interpolY1 = gtetaFunctionDerived[itype][interpolIDX];
        interpolY2 = gtetaFunctionDerived[itype][interpolIDX+1];
        gtetaFunctionDerived_temp = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

        preGtetaFunction[neighbor_j][neighbor_k]=params[ijkparam].gamma*gtetaFunction_temp;
        preGtetaFunctionDerived[neighbor_j][neighbor_k]=params[ijkparam].gamma*gtetaFunctionDerived_temp;
        preGtetaFunction[neighbor_k][neighbor_j]=params[ijkparam].gamma*gtetaFunction_temp;
        preGtetaFunctionDerived[neighbor_k][neighbor_j]=params[ijkparam].gamma*gtetaFunctionDerived_temp;

      } // loop on K

    } // loop on J


    // loop over neighbors of atom i
    for (int neighbor_j = 0; neighbor_j < jnum; neighbor_j++) {

      double dr_ij[3], r_ij, f_ij[3];

      j = jlist[neighbor_j];
      j &= NEIGHMASK;

      dr_ij[0] = xtmp - x[j][0];
      dr_ij[1] = ytmp - x[j][1];
      dr_ij[2] = ztmp - x[j][2];
      r_ij = dr_ij[0]*dr_ij[0] + dr_ij[1]*dr_ij[1] + dr_ij[2]*dr_ij[2];

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];

      if (r_ij > params[ijparam].cutsq) continue;

      r_ij = sqrt(r_ij);
      invR_ij = 1.0 / r_ij;

      directorCos_ij_x = invR_ij * dr_ij[0];
      directorCos_ij_y = invR_ij * dr_ij[1];
      directorCos_ij_z = invR_ij * dr_ij[2];

      exponentRepulsivePotential = params[ijparam].lam1 * r_ij;
      exponentAttractivePotential = params[ijparam].lam2 * r_ij;

      // repulsiveExponential
      interpolDeltaX =  exponentRepulsivePotential - minArgumentExponential;
      interpolTMP = (interpolDeltaX * GRIDDENSITY_EXP);
      interpolIDX = (int) interpolTMP;
      interpolY1 = exponential[interpolIDX];
      interpolY2 = exponential[interpolIDX+1];
      repulsiveExponential = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);
      // attractiveExponential
      interpolDeltaX =  exponentAttractivePotential - minArgumentExponential;
      interpolTMP = (interpolDeltaX * GRIDDENSITY_EXP);
      interpolIDX = (int) interpolTMP;
      interpolY1 = exponential[interpolIDX];
      interpolY2 = exponential[interpolIDX+1];
      attractiveExponential = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

      repulsivePotential = params[ijparam].biga * repulsiveExponential;
      attractivePotential = -params[ijparam].bigb * attractiveExponential;

      cutoffFunctionIJ = preCutoffFunction[neighbor_j];
      cutoffFunctionDerivedIJ = preCutoffFunctionDerived[neighbor_j];

      zeta = 0.0;

      // first loop over neighbors of atom i except j - part 1/2
      for (int neighbor_k = 0; neighbor_k < neighbor_j; neighbor_k++) {
        double dr_ik[3], r_ik;

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = xtmp -x[k][0];
        dr_ik[1] = ytmp -x[k][1];
        dr_ik[2] = ztmp -x[k][2];
        r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        gtetaFunctionIJK = preGtetaFunction[neighbor_j][neighbor_k];

        cutoffFunctionIK = preCutoffFunction[neighbor_k];

        zeta += cutoffFunctionIK * gtetaFunctionIJK;

      }

      // first loop over neighbors of atom i except j - part 2/2
      for (int neighbor_k = neighbor_j+1; neighbor_k < jnum; neighbor_k++) {
        double dr_ik[3], r_ik;

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = xtmp -x[k][0];
        dr_ik[1] = ytmp -x[k][1];
        dr_ik[2] = ztmp -x[k][2];
        r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        gtetaFunctionIJK = preGtetaFunction[neighbor_j][neighbor_k];

        cutoffFunctionIK = preCutoffFunction[neighbor_k];

        zeta += cutoffFunctionIK * gtetaFunctionIJK;
      }

      // betaZetaPowerIJK
      interpolDeltaX= params[ijparam].beta * zeta;
      interpolTMP = (interpolDeltaX * GRIDDENSITY_BIJ);
      interpolIDX = (int) interpolTMP;
      interpolY1 = betaZetaPower[itype][interpolIDX];
      interpolY2 = betaZetaPower[itype][interpolIDX+1];
      betaZetaPowerIJK = (interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX));
      // betaZetaPowerDerivedIJK
      interpolY1 = betaZetaPowerDerived[itype][interpolIDX];
      interpolY2 = betaZetaPowerDerived[itype][interpolIDX+1];
      betaZetaPowerDerivedIJK = params[ijparam].beta*(interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX));

      // Forces and virial
      factor_force_ij = 0.5*cutoffFunctionDerivedIJ*(repulsivePotential + attractivePotential * betaZetaPowerIJK)+0.5*cutoffFunctionIJ*(-repulsivePotential*params[ijparam].lam1-betaZetaPowerIJK*attractivePotential*params[ijparam].lam2);

      f_ij[0] = factor_force_ij * directorCos_ij_x;
      f_ij[1] = factor_force_ij * directorCos_ij_y;
      f_ij[2] = factor_force_ij * directorCos_ij_z;

      f[j][0] += f_ij[0];
      f[j][1] += f_ij[1];
      f[j][2] += f_ij[2];

      fxtmp -= f_ij[0];
      fytmp -= f_ij[1];
      fztmp -= f_ij[2];

      // potential energy
      evdwl = cutoffFunctionIJ * repulsivePotential + cutoffFunctionIJ * attractivePotential * betaZetaPowerIJK;

      if (evflag) ev_tally(i, j, nlocal, newton_pair, 0.5 * evdwl, 0.0,
                           -factor_force_ij*invR_ij, dr_ij[0], dr_ij[1], dr_ij[2]);

      factor_force_tot= 0.5*cutoffFunctionIJ*attractivePotential*betaZetaPowerDerivedIJK;

      // second loop over neighbors of atom i except j, forces and virial only - part 1/2
      for (int neighbor_k = 0; neighbor_k < neighbor_j; neighbor_k++) {
        double dr_ik[3], r_ik, f_ik[3];

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = xtmp -x[k][0];
        dr_ik[1] = ytmp -x[k][1];
        dr_ik[2] = ztmp -x[k][2];
        r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        r_ik = sqrt(r_ik);
        invR_ik = 1.0 / r_ik;

        directorCos_ik_x = invR_ik * dr_ik[0];
        directorCos_ik_y = invR_ik * dr_ik[1];
        directorCos_ik_z = invR_ik * dr_ik[2];

        cosTeta = directorCos_ij_x * directorCos_ik_x + directorCos_ij_y * directorCos_ik_y + directorCos_ij_z * directorCos_ik_z;

        gtetaFunctionIJK = preGtetaFunction[neighbor_j][neighbor_k];

        gtetaFunctionDerivedIJK = preGtetaFunctionDerived[neighbor_j][neighbor_k];

        cutoffFunctionIK = preCutoffFunction[neighbor_k];

        cutoffFunctionDerivedIK = preCutoffFunctionDerived[neighbor_k];

        factor_force3_ij= cutoffFunctionIK * gtetaFunctionDerivedIJK * invR_ij *factor_force_tot;

        f_ij[0] = factor_force3_ij * (directorCos_ij_x*cosTeta - directorCos_ik_x);
        f_ij[1] = factor_force3_ij * (directorCos_ij_y*cosTeta - directorCos_ik_y);
        f_ij[2] = factor_force3_ij * (directorCos_ij_z*cosTeta - directorCos_ik_z);

        factor_1_force3_ik = (cutoffFunctionIK * gtetaFunctionDerivedIJK * invR_ik)*factor_force_tot;
        factor_2_force3_ik = -(cutoffFunctionDerivedIK * gtetaFunctionIJK)*factor_force_tot;

        f_ik[0] = factor_1_force3_ik * (directorCos_ik_x*cosTeta - directorCos_ij_x) + factor_2_force3_ik * directorCos_ik_x;
        f_ik[1] = factor_1_force3_ik * (directorCos_ik_y*cosTeta - directorCos_ij_y) + factor_2_force3_ik * directorCos_ik_y;
        f_ik[2] = factor_1_force3_ik * (directorCos_ik_z*cosTeta - directorCos_ij_z) + factor_2_force3_ik * directorCos_ik_z;

        f[j][0] -= f_ij[0];
        f[j][1] -= f_ij[1];
        f[j][2] -= f_ij[2];

        f[k][0] -= f_ik[0];
        f[k][1] -= f_ik[1];
        f[k][2] -= f_ik[2];

        fxtmp += f_ij[0] + f_ik[0];
        fytmp += f_ij[1] + f_ik[1];
        fztmp += f_ij[2] + f_ik[2];

        if (vflag_either) v_tally3(i,j,k,f_ij,f_ik,dr_ij,dr_ik);
      }

      // second loop over neighbors of atom i except j, forces and virial only - part 2/2
      for (int neighbor_k = neighbor_j+1; neighbor_k < jnum; neighbor_k++) {
        double dr_ik[3], r_ik, f_ik[3];

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = xtmp -x[k][0];
        dr_ik[1] = ytmp -x[k][1];
        dr_ik[2] = ztmp -x[k][2];
        r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        r_ik = sqrt(r_ik);
        invR_ik = 1.0 / r_ik;

        directorCos_ik_x = invR_ik * dr_ik[0];
        directorCos_ik_y = invR_ik * dr_ik[1];
        directorCos_ik_z = invR_ik * dr_ik[2];

        cosTeta = directorCos_ij_x * directorCos_ik_x + directorCos_ij_y * directorCos_ik_y + directorCos_ij_z * directorCos_ik_z;

        gtetaFunctionIJK = preGtetaFunction[neighbor_j][neighbor_k];

        gtetaFunctionDerivedIJK = preGtetaFunctionDerived[neighbor_j][neighbor_k];

        cutoffFunctionIK = preCutoffFunction[neighbor_k];

        cutoffFunctionDerivedIK = preCutoffFunctionDerived[neighbor_k];

        factor_force3_ij= cutoffFunctionIK * gtetaFunctionDerivedIJK * invR_ij *factor_force_tot;

        f_ij[0] = factor_force3_ij * (directorCos_ij_x*cosTeta - directorCos_ik_x);
        f_ij[1] = factor_force3_ij * (directorCos_ij_y*cosTeta - directorCos_ik_y);
        f_ij[2] = factor_force3_ij * (directorCos_ij_z*cosTeta - directorCos_ik_z);

        factor_1_force3_ik = (cutoffFunctionIK * gtetaFunctionDerivedIJK * invR_ik)*factor_force_tot;
        factor_2_force3_ik = -(cutoffFunctionDerivedIK * gtetaFunctionIJK)*factor_force_tot;

        f_ik[0] = factor_1_force3_ik * (directorCos_ik_x*cosTeta - directorCos_ij_x) + factor_2_force3_ik * directorCos_ik_x;
        f_ik[1] = factor_1_force3_ik * (directorCos_ik_y*cosTeta - directorCos_ij_y) + factor_2_force3_ik * directorCos_ik_y;
        f_ik[2] = factor_1_force3_ik * (directorCos_ik_z*cosTeta - directorCos_ij_z) + factor_2_force3_ik * directorCos_ik_z;

        f[j][0] -= f_ij[0];
        f[j][1] -= f_ij[1];
        f[j][2] -= f_ij[2];

        f[k][0] -= f_ik[0];
        f[k][1] -= f_ik[1];
        f[k][2] -= f_ik[2];

        fxtmp += f_ij[0] + f_ik[0];
        fytmp += f_ij[1] + f_ik[1];
        fztmp += f_ij[2] + f_ik[2];

        if (vflag_either) v_tally3(i,j,k,f_ij,f_ik,dr_ij,dr_ik);
      }
    } // loop on J
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  } // loop on I

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairTersoffTable::deallocatePreLoops(void)
{
    memory->destroy(preGtetaFunction);
    memory->destroy(preGtetaFunctionDerived);
    memory->destroy(preCutoffFunction);
    memory->destroy(preCutoffFunctionDerived);
}

void PairTersoffTable::allocatePreLoops(void)
{
  deallocatePreLoops();
  memory->create(preGtetaFunction,leadingDimensionInteractionList,
                 leadingDimensionInteractionList,"tersofftable:preGtetaFunction");
  memory->create(preGtetaFunctionDerived,leadingDimensionInteractionList,
                 leadingDimensionInteractionList,"tersofftable:preGtetaFunctionDerived");
  memory->create(preCutoffFunction,leadingDimensionInteractionList,
                 "tersofftable:preCutoffFunction");
  memory->create(preCutoffFunctionDerived,leadingDimensionInteractionList,
                 "tersofftable:preCutoffFunctionDerived");
}

void PairTersoffTable::deallocateGrids()
{
  memory->destroy(exponential);
  memory->destroy(gtetaFunction);
  memory->destroy(gtetaFunctionDerived);
  memory->destroy(cutoffFunction);
  memory->destroy(cutoffFunctionDerived);
  memory->destroy(betaZetaPower);
  memory->destroy(betaZetaPowerDerived);
}

void PairTersoffTable::allocateGrids(void)
{
  int   i, j, k, l;

  int     numGridPointsExponential, numGridPointsGtetaFunction, numGridPointsOneCutoffFunction;
  int     numGridPointsNotOneCutoffFunction, numGridPointsCutoffFunction, numGridPointsBetaZetaPower;
  // double minArgumentExponential;
  double  deltaArgumentCutoffFunction, deltaArgumentExponential, deltaArgumentBetaZetaPower;
  double  deltaArgumentGtetaFunction;
  double  r, minMu, maxLambda, maxCutoff;
  double const PI=acos(-1.0);

  deallocateGrids();

  // exponential

  // find min and max argument
  minMu=params[0].lam2;
  maxLambda=params[0].lam1;
  for (i=1; i<nparams; i++) {
    if (params[i].lam2 < minMu) minMu = params[i].lam2;
    if (params[i].lam1 > maxLambda) maxLambda = params[i].lam1;
  }
  maxCutoff=cutmax;

  minArgumentExponential=minMu*GRIDSTART;
  numGridPointsExponential=(int)((maxLambda*maxCutoff-minArgumentExponential)*GRIDDENSITY_EXP)+2;
  memory->create(exponential,numGridPointsExponential,"tersofftable:exponential");

  r = minArgumentExponential;
  deltaArgumentExponential = 1.0 / GRIDDENSITY_EXP;
  for (i = 0; i < numGridPointsExponential; i++)
    {
      exponential[i] = exp(-r);
      r += deltaArgumentExponential;
    }


  // gtetaFunction

  numGridPointsGtetaFunction=(int)(2.0*GRIDDENSITY_GTETA)+2;

  memory->create(gtetaFunction,nelements,numGridPointsGtetaFunction,"tersofftable:gtetaFunction");
  memory->create(gtetaFunctionDerived,nelements,numGridPointsGtetaFunction,"tersofftable:gtetaFunctionDerived");

  r = minArgumentExponential;
  for (i=0; i<nelements; i++) {
    r = -1.0;
    deltaArgumentGtetaFunction = 1.0 / GRIDDENSITY_GTETA;

    int iparam = elem3param[i][i][i];
    double c = params[iparam].c;
    double d = params[iparam].d;
    double h = params[iparam].h;

    for (j = 0; j < numGridPointsGtetaFunction; j++) {
      gtetaFunction[i][j]=1.0+(c*c)/(d*d)-(c*c)/(d*d+(h-r)*(h-r));
      gtetaFunctionDerived[i][j]= -2.0 * c * c * (h-r) / ((d*d+(h-r)*(h-r))*(d*d+(h-r)*(h-r)));
      r += deltaArgumentGtetaFunction;
    }
  }


  // cutoffFunction, zetaFunction, find grids.

  int ngrid_max = -1;
  int zeta_max = -1;

  for (i=0; i<nelements; i++) {

    int iparam = elem3param[i][i][i];
    double c = params[iparam].c;
    double d = params[iparam].d;
    double beta = params[iparam].beta;

    numGridPointsBetaZetaPower=(int)(((1.0+(c*c)/(d*d)-(c*c)/(d*d+4))*beta*leadingDimensionInteractionList*GRIDDENSITY_BIJ))+2;
    zeta_max = MAX(zeta_max,numGridPointsBetaZetaPower);

    for (j=0; j<nelements; j++) {
      for (k=0; k<nelements; k++) {

        int ijparam = elem3param[i][j][k];
        double cutoffR = params[ijparam].cutoffR;
        double cutoffS = params[ijparam].cutoffS;

        numGridPointsOneCutoffFunction=(int) ((cutoffR-GRIDSTART)*GRIDDENSITY_FCUTOFF)+1;
        numGridPointsNotOneCutoffFunction=(int) ((cutoffS-cutoffR)*GRIDDENSITY_FCUTOFF)+2;
        numGridPointsCutoffFunction=numGridPointsOneCutoffFunction+numGridPointsNotOneCutoffFunction;

        ngrid_max = MAX(ngrid_max,numGridPointsCutoffFunction);
      }
    }
  }

  memory->create(cutoffFunction,nelements,nelements,ngrid_max,"tersoff:cutfunc");
  memory->create(cutoffFunctionDerived,nelements,nelements,ngrid_max,"tersoff:cutfuncD");

  // cutoffFunction, compute.

  for (i=0; i<nelements; i++) {
    for (j=0; j<nelements; j++) {
      for (j=0; j<nelements; j++) {
        int ijparam = elem3param[i][j][j];
        double cutoffR = params[ijparam].cutoffR;
        double cutoffS = params[ijparam].cutoffS;

        numGridPointsOneCutoffFunction=(int) ((cutoffR-GRIDSTART)*GRIDDENSITY_FCUTOFF)+1;
        numGridPointsNotOneCutoffFunction=(int) ((cutoffS-cutoffR)*GRIDDENSITY_FCUTOFF)+2;
        numGridPointsCutoffFunction=numGridPointsOneCutoffFunction+numGridPointsNotOneCutoffFunction;

        r = GRIDSTART;
        deltaArgumentCutoffFunction = 1.0 / GRIDDENSITY_FCUTOFF;

        for (l = 0; l < numGridPointsOneCutoffFunction; l++) {
          cutoffFunction[i][j][l] = 1.0;
          cutoffFunctionDerived[i][j][l]=0.0;
          r += deltaArgumentCutoffFunction;
        }

        for (l = numGridPointsOneCutoffFunction; l < numGridPointsCutoffFunction; l++) {
          cutoffFunction[i][j][l] = 0.5 + 0.5 * cos (PI * (r - cutoffR)/(cutoffS-cutoffR)) ;
          cutoffFunctionDerived[i][j][l] =  -0.5 * PI * sin (PI * (r - cutoffR)/(cutoffS-cutoffR)) / (cutoffS-cutoffR) ;
          r += deltaArgumentCutoffFunction;
        }
      }
    }
  }

  // betaZetaPower, compute

  memory->create(betaZetaPower,nelements,zeta_max,"tersoff:zetafunc");
  memory->create(betaZetaPowerDerived,nelements,zeta_max,"tersoff:zetafuncD");

  for (i=0; i<nelements; i++) {

    int iparam = elem3param[i][i][i];
    double c = params[iparam].c;
    double d = params[iparam].d;
    double beta = params[iparam].beta;

    numGridPointsBetaZetaPower=(int)(((1.0+(c*c)/(d*d)-(c*c)/(d*d+4))*beta*leadingDimensionInteractionList*GRIDDENSITY_BIJ))+2;

    r=0.0;
    deltaArgumentBetaZetaPower = 1.0 / GRIDDENSITY_BIJ;

    betaZetaPower[i][0]=1.0;

    r += deltaArgumentBetaZetaPower;

    for (j = 1; j < numGridPointsBetaZetaPower; j++) {
      double powern=params[iparam].powern;
      betaZetaPower[i][j]=pow((1+pow(r,powern)),-1/(2*powern));
      betaZetaPowerDerived[i][j]=-0.5*pow(r,powern-1.0)*pow((1+pow(r,powern)),-1/(2*powern)-1) ;
      r += deltaArgumentBetaZetaPower;
    }
    betaZetaPowerDerived[i][0]=(betaZetaPower[i][1]-1.0)*GRIDDENSITY_BIJ;
  }
}

void PairTersoffTable::allocate()
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

void PairTersoffTable::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTersoffTable::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg-3,arg+3);

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();

  // allocate tables and internal structures

  allocatePreLoops();
  allocateGrids();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTersoffTable::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style tersoff/table requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTersoffTable::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairTersoffTable::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "tersoff/table", unit_convert_flag);
    char *line;

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY,
                                                            unit_convert);

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

        // some parameters are not used since only Tersoff_2 is implemented

        params[nparams].ielement  = ielement;
        params[nparams].jelement  = jelement;
        params[nparams].kelement  = kelement;
        params[nparams].powerm    = values.next_double(); // not used
        params[nparams].gamma     = values.next_double(); // not used
        params[nparams].lam3      = values.next_double(); // not used
        params[nparams].c         = values.next_double();
        params[nparams].d         = values.next_double();
        params[nparams].h         = values.next_double();
        params[nparams].powern    = values.next_double();
        params[nparams].beta      = values.next_double();
        params[nparams].lam2      = values.next_double();
        params[nparams].bigb      = values.next_double();
        double bigr               = values.next_double();
        double bigd               = values.next_double();
        params[nparams].cutoffR = bigr - bigd;
        params[nparams].cutoffS = bigr + bigd;
        params[nparams].lam1      = values.next_double();
        params[nparams].biga      = values.next_double();

        if (unit_convert) {
          params[nparams].biga *= conversion_factor;
          params[nparams].bigb *= conversion_factor;
        }
      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      if (params[nparams].c < 0.0 ||
          params[nparams].d < 0.0 ||
          params[nparams].powern < 0.0 ||
          params[nparams].beta < 0.0 ||
          params[nparams].lam2 < 0.0 ||
          params[nparams].bigb < 0.0 ||
          params[nparams].cutoffR < 0.0 ||
          params[nparams].cutoffS < 0.0 ||
          params[nparams].cutoffR > params[nparams].cutoffS ||
          params[nparams].lam1 < 0.0 ||
          params[nparams].biga < 0.0
        ) error->one(FLERR,"Illegal Tersoff parameter");

      // only tersoff_2 parametrization is implemented

      if (params[nparams].gamma != 1.0 || params[nparams].lam3 != 0.0)
        error->one(FLERR,"Currently the tersoff/table pair_style only "
                   "implements the Tersoff_2 parametrization");
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

void PairTersoffTable::setup_params()
{
  int i,j,k,m,n;

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
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem3param[i][j][k] = n;
      }

  // set cutoff square
  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].cutoffS;
    params[m].cutsq = params[m].cut*params[m].cut;
  }

  // set cutmax to max of all params
  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    if (params[m].cut > cutmax) cutmax = params[m].cut;
  }
}

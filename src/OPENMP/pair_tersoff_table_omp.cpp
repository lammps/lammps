// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_tersoff_table_omp.h"

#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

#define GRIDSTART 0.1
#define GRIDDENSITY_FCUTOFF 5000
#define GRIDDENSITY_EXP 12000
#define GRIDDENSITY_GTETA 12000
#define GRIDDENSITY_BIJ 7500

// max number of interaction per atom for environment potential

#define leadingDimensionInteractionList 64

/* ---------------------------------------------------------------------- */

PairTersoffTableOMP::PairTersoffTableOMP(LAMMPS *lmp) :
  PairTersoffTable(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;

  thrGtetaFunction = thrGtetaFunctionDerived = nullptr;
  thrCutoffFunction = thrCutoffFunctionDerived = nullptr;
}

/* ---------------------------------------------------------------------- */

PairTersoffTableOMP::~PairTersoffTableOMP()
{
  if (allocated) {
    PairTersoffTableOMP::deallocatePreLoops();
  }
}


/* ---------------------------------------------------------------------- */

void PairTersoffTableOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag)
      if (vflag_either) eval<1,1>(ifrom, ito, thr);
      else eval<1,0>(ifrom, ito, thr);
    else eval<0,0>(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int VFLAG_EITHER>
void PairTersoffTableOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,k,ii,jnum;
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

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const int tid = thr->get_tid();

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms
  for (ii = iifrom; ii < iito; ii++) {

    i = ilist[ii];
    itype = map[type[i]];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    fxtmp = fytmp = fztmp = 0.0;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (check_error_thr((jnum > leadingDimensionInteractionList), tid,
                        FLERR,"Too many neighbors for interaction list.\n"
                        "Check your system or increase 'leadingDimension"
                        "InteractionList'"))
      return;

    // Pre-calculate gteta and cutoff function
    for (int neighbor_j = 0; neighbor_j < jnum; neighbor_j++) {

      double dr_ij[3], r_ij;

      j = jlist[neighbor_j];
      j &= NEIGHMASK;

      dr_ij[0] = xtmp - x[j].x;
      dr_ij[1] = ytmp - x[j].y;
      dr_ij[2] = ztmp - x[j].z;
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
      thrCutoffFunction[tid][neighbor_j] = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);
      // preCutoffFunctionDerived
      interpolY1 = cutoffFunctionDerived[itype][jtype][interpolIDX];
      interpolY2 = cutoffFunctionDerived[itype][jtype][interpolIDX+1];
      thrCutoffFunctionDerived[tid][neighbor_j] = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);


      for (int neighbor_k = neighbor_j + 1; neighbor_k < jnum; neighbor_k++) {
        double dr_ik[3], r_ik;

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = xtmp -x[k].x;
        dr_ik[1] = ytmp -x[k].y;
        dr_ik[2] = ztmp -x[k].z;
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

        thrGtetaFunction[tid][neighbor_j][neighbor_k]=params[ijkparam].gamma*gtetaFunction_temp;
        thrGtetaFunctionDerived[tid][neighbor_j][neighbor_k]=params[ijkparam].gamma*gtetaFunctionDerived_temp;
        thrGtetaFunction[tid][neighbor_k][neighbor_j]=params[ijkparam].gamma*gtetaFunction_temp;
        thrGtetaFunctionDerived[tid][neighbor_k][neighbor_j]=params[ijkparam].gamma*gtetaFunctionDerived_temp;

      } // loop on K

    } // loop on J


    // loop over neighbors of atom i
    for (int neighbor_j = 0; neighbor_j < jnum; neighbor_j++) {

      double dr_ij[3], r_ij, f_ij[3];

      j = jlist[neighbor_j];
      j &= NEIGHMASK;

      dr_ij[0] = xtmp - x[j].x;
      dr_ij[1] = ytmp - x[j].y;
      dr_ij[2] = ztmp - x[j].z;
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

      cutoffFunctionIJ = thrCutoffFunction[tid][neighbor_j];
      cutoffFunctionDerivedIJ = thrCutoffFunctionDerived[tid][neighbor_j];

      zeta = 0.0;

      // first loop over neighbors of atom i except j - part 1/2
      for (int neighbor_k = 0; neighbor_k < neighbor_j; neighbor_k++) {
        double dr_ik[3], r_ik;

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = xtmp -x[k].x;
        dr_ik[1] = ytmp -x[k].y;
        dr_ik[2] = ztmp -x[k].z;
        r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        gtetaFunctionIJK = thrGtetaFunction[tid][neighbor_j][neighbor_k];

        cutoffFunctionIK = thrCutoffFunction[tid][neighbor_k];

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

        dr_ik[0] = xtmp -x[k].x;
        dr_ik[1] = ytmp -x[k].y;
        dr_ik[2] = ztmp -x[k].z;
        r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        gtetaFunctionIJK = thrGtetaFunction[tid][neighbor_j][neighbor_k];

        cutoffFunctionIK = thrCutoffFunction[tid][neighbor_k];

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

      f[j].x += f_ij[0];
      f[j].y += f_ij[1];
      f[j].z += f_ij[2];

      fxtmp -= f_ij[0];
      fytmp -= f_ij[1];
      fztmp -= f_ij[2];

      // potential energy
      evdwl = cutoffFunctionIJ * repulsivePotential
        + cutoffFunctionIJ * attractivePotential * betaZetaPowerIJK;

      if (EVFLAG) ev_tally_thr(this,i, j, nlocal, /* newton_pair */ 1, 0.5 * evdwl, 0.0,
                               -factor_force_ij*invR_ij, dr_ij[0], dr_ij[1], dr_ij[2],thr);

      factor_force_tot= 0.5*cutoffFunctionIJ*attractivePotential*betaZetaPowerDerivedIJK;

      // second loop over neighbors of atom i except j, forces and virial only - part 1/2
      for (int neighbor_k = 0; neighbor_k < neighbor_j; neighbor_k++) {
        double dr_ik[3], r_ik, f_ik[3];

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = xtmp -x[k].x;
        dr_ik[1] = ytmp -x[k].y;
        dr_ik[2] = ztmp -x[k].z;
        r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        r_ik = sqrt(r_ik);
        invR_ik = 1.0 / r_ik;

        directorCos_ik_x = invR_ik * dr_ik[0];
        directorCos_ik_y = invR_ik * dr_ik[1];
        directorCos_ik_z = invR_ik * dr_ik[2];

        cosTeta = directorCos_ij_x * directorCos_ik_x + directorCos_ij_y * directorCos_ik_y
          + directorCos_ij_z * directorCos_ik_z;

        gtetaFunctionIJK = thrGtetaFunction[tid][neighbor_j][neighbor_k];

        gtetaFunctionDerivedIJK = thrGtetaFunctionDerived[tid][neighbor_j][neighbor_k];

        cutoffFunctionIK = thrCutoffFunction[tid][neighbor_k];

        cutoffFunctionDerivedIK = thrCutoffFunctionDerived[tid][neighbor_k];

        factor_force3_ij= cutoffFunctionIK * gtetaFunctionDerivedIJK * invR_ij *factor_force_tot;

        f_ij[0] = factor_force3_ij * (directorCos_ij_x*cosTeta - directorCos_ik_x);
        f_ij[1] = factor_force3_ij * (directorCos_ij_y*cosTeta - directorCos_ik_y);
        f_ij[2] = factor_force3_ij * (directorCos_ij_z*cosTeta - directorCos_ik_z);

        factor_1_force3_ik = (cutoffFunctionIK * gtetaFunctionDerivedIJK * invR_ik)*factor_force_tot;
        factor_2_force3_ik = -(cutoffFunctionDerivedIK * gtetaFunctionIJK)*factor_force_tot;

        f_ik[0] = factor_1_force3_ik * (directorCos_ik_x*cosTeta - directorCos_ij_x)
          + factor_2_force3_ik * directorCos_ik_x;
        f_ik[1] = factor_1_force3_ik * (directorCos_ik_y*cosTeta - directorCos_ij_y)
          + factor_2_force3_ik * directorCos_ik_y;
        f_ik[2] = factor_1_force3_ik * (directorCos_ik_z*cosTeta - directorCos_ij_z)
          + factor_2_force3_ik * directorCos_ik_z;

        f[j].x -= f_ij[0];
        f[j].y -= f_ij[1];
        f[j].z -= f_ij[2];

        f[k].x -= f_ik[0];
        f[k].y -= f_ik[1];
        f[k].z -= f_ik[2];

        fxtmp += f_ij[0] + f_ik[0];
        fytmp += f_ij[1] + f_ik[1];
        fztmp += f_ij[2] + f_ik[2];

        if (VFLAG_EITHER) v_tally3_thr(this,i,j,k,f_ij,f_ik,dr_ij,dr_ik,thr);
      }

      // second loop over neighbors of atom i except j, forces and virial only - part 2/2
      for (int neighbor_k = neighbor_j+1; neighbor_k < jnum; neighbor_k++) {
        double dr_ik[3], r_ik, f_ik[3];

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = xtmp -x[k].x;
        dr_ik[1] = ytmp -x[k].y;
        dr_ik[2] = ztmp -x[k].z;
        r_ik = dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        r_ik = sqrt(r_ik);
        invR_ik = 1.0 / r_ik;

        directorCos_ik_x = invR_ik * dr_ik[0];
        directorCos_ik_y = invR_ik * dr_ik[1];
        directorCos_ik_z = invR_ik * dr_ik[2];

        cosTeta = directorCos_ij_x * directorCos_ik_x + directorCos_ij_y * directorCos_ik_y
          + directorCos_ij_z * directorCos_ik_z;

        gtetaFunctionIJK = thrGtetaFunction[tid][neighbor_j][neighbor_k];

        gtetaFunctionDerivedIJK = thrGtetaFunctionDerived[tid][neighbor_j][neighbor_k];

        cutoffFunctionIK = thrCutoffFunction[tid][neighbor_k];

        cutoffFunctionDerivedIK = thrCutoffFunctionDerived[tid][neighbor_k];

        factor_force3_ij= cutoffFunctionIK * gtetaFunctionDerivedIJK * invR_ij *factor_force_tot;

        f_ij[0] = factor_force3_ij * (directorCos_ij_x*cosTeta - directorCos_ik_x);
        f_ij[1] = factor_force3_ij * (directorCos_ij_y*cosTeta - directorCos_ik_y);
        f_ij[2] = factor_force3_ij * (directorCos_ij_z*cosTeta - directorCos_ik_z);

        factor_1_force3_ik = (cutoffFunctionIK * gtetaFunctionDerivedIJK * invR_ik)*factor_force_tot;
        factor_2_force3_ik = -(cutoffFunctionDerivedIK * gtetaFunctionIJK)*factor_force_tot;

        f_ik[0] = factor_1_force3_ik * (directorCos_ik_x*cosTeta - directorCos_ij_x)
          + factor_2_force3_ik * directorCos_ik_x;
        f_ik[1] = factor_1_force3_ik * (directorCos_ik_y*cosTeta - directorCos_ij_y)
          + factor_2_force3_ik * directorCos_ik_y;
        f_ik[2] = factor_1_force3_ik * (directorCos_ik_z*cosTeta - directorCos_ij_z)
          + factor_2_force3_ik * directorCos_ik_z;

        f[j].x -= f_ij[0];
        f[j].y -= f_ij[1];
        f[j].z -= f_ij[2];

        f[k].x -= f_ik[0];
        f[k].y -= f_ik[1];
        f[k].z -= f_ik[2];

        fxtmp += f_ij[0] + f_ik[0];
        fytmp += f_ij[1] + f_ik[1];
        fztmp += f_ij[2] + f_ik[2];

        if (VFLAG_EITHER) v_tally3_thr(this,i,j,k,f_ij,f_ik,dr_ij,dr_ik,thr);
      }
    } // loop on J
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  } // loop on I
}

void PairTersoffTableOMP::deallocatePreLoops()
{
    memory->destroy(thrGtetaFunction);
    memory->destroy(thrGtetaFunctionDerived);
    memory->destroy(thrCutoffFunction);
    memory->destroy(thrCutoffFunctionDerived);
}

void PairTersoffTableOMP::allocatePreLoops()
{
  const int nthreads = comm->nthreads;
  memory->create(thrGtetaFunction,nthreads,leadingDimensionInteractionList,leadingDimensionInteractionList,"tersofftable:thrGtetaFunction");

  memory->create(thrGtetaFunctionDerived,nthreads,leadingDimensionInteractionList,leadingDimensionInteractionList,"tersofftable:thrGtetaFunctionDerived");

  memory->create(thrCutoffFunction,nthreads,leadingDimensionInteractionList,"tersofftable:thrCutoffFunction");

  memory->create(thrCutoffFunctionDerived,nthreads,leadingDimensionInteractionList,"tersofftable:thrCutoffFunctionDerived");
}

/* ---------------------------------------------------------------------- */

double PairTersoffTableOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairTersoffTable::memory_usage();

  return bytes;
}

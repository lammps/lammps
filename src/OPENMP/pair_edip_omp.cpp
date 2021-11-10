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

#include "pair_edip_omp.h"

#include "atom.h"
#include "comm.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

// max number of interaction per atom for f(Z) environment potential

static constexpr int leadDimInteractionList = 64;

#define GRIDDENSITY 8000
#define GRIDSTART 0.1

/* ---------------------------------------------------------------------- */

PairEDIPOMP::PairEDIPOMP(LAMMPS *lmp) : PairEDIP(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairEDIPOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag) {
      if (eflag) {
        if (vflag_atom)
          eval<1, 1, 1>(ifrom, ito, thr);
        else
          eval<1, 1, 0>(ifrom, ito, thr);
      } else {
        if (vflag_atom)
          eval<1, 0, 1>(ifrom, ito, thr);
        else
          eval<1, 0, 0>(ifrom, ito, thr);
      }
    } else
      eval<0, 0, 0>(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int VFLAG_ATOM>
void PairEDIPOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  int i, j, k, ii, jnum;
  int itype, jtype, ktype, ijparam, ikparam;
  double xtmp, ytmp, ztmp, evdwl;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int preForceCoord_counter;

  double invR_ij;
  double invR_ik;
  double directorCos_ij_x;
  double directorCos_ij_y;
  double directorCos_ij_z;
  double directorCos_ik_x;
  double directorCos_ik_y;
  double directorCos_ik_z;
  double cosTeta;

  int interpolIDX;
  double interpolTMP;
  double interpolDeltaX;
  double interpolY1;
  double interpolY2;

  double invRMinusCutoffA;
  double sigmaInvRMinusCutoffA;
  double gammInvRMinusCutoffA;
  double cosTetaDiff;
  double cosTetaDiffCosTetaDiff;
  double cutoffFunction_ij;
  double exp2B_ij;
  double exp2BDerived_ij;
  double pow2B_ij;
  double pow2BDerived_ij;
  double exp3B_ij;
  double exp3BDerived_ij;
  double exp3B_ik;
  double exp3BDerived_ik;
  double qFunction;
  double tauFunction;
  double tauFunctionDerived;
  double expMinusBetaZeta_iZeta_i;
  double qFunctionCosTetaDiffCosTetaDiff;
  double expMinusQFunctionCosTetaDiffCosTetaDiff;
  double zeta_i;
  double zeta_iDerived;
  double zeta_iDerivedInvR_ij;

  double forceModCoord_factor;
  double forceModCoord;
  double forceModCoord_ij;
  double forceMod2B;
  double forceMod3B_factor1_ij;
  double forceMod3B_factor2_ij;
  double forceMod3B_factor2;
  double forceMod3B_factor1_ik;
  double forceMod3B_factor2_ik;
  double potentia3B_factor;
  double potential2B_factor;

  const int tid = thr->get_tid();

  double *pre_thrInvR_ij = preInvR_ij + tid * leadDimInteractionList;
  double *pre_thrExp3B_ij = preExp3B_ij + tid * leadDimInteractionList;
  double *pre_thrExp3BDerived_ij = preExp3BDerived_ij + tid * leadDimInteractionList;
  double *pre_thrExp2B_ij = preExp2B_ij + tid * leadDimInteractionList;
  double *pre_thrExp2BDerived_ij = preExp2BDerived_ij + tid * leadDimInteractionList;
  double *pre_thrPow2B_ij = prePow2B_ij + tid * leadDimInteractionList;
  double *pre_thrForceCoord = preForceCoord + tid * leadDimInteractionList;

  const dbl3_t *_noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int *_noalias const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms

  for (ii = iifrom; ii < iito; ii++) {
    zeta_i = 0.0;
    int numForceCoordPairs = 0;

    i = ilist[ii];
    itype = map[type[i]];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // pre-loop to compute environment coordination f(Z)

    for (int neighbor_j = 0; neighbor_j < jnum; neighbor_j++) {
      j = jlist[neighbor_j];
      j &= NEIGHMASK;

      double dr_ij[3], r_ij;

      dr_ij[0] = xtmp - x[j].x;
      dr_ij[1] = ytmp - x[j].y;
      dr_ij[2] = ztmp - x[j].z;
      r_ij = dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1] + dr_ij[2] * dr_ij[2];

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (r_ij > params[ijparam].cutsq) continue;

      r_ij = sqrt(r_ij);

      invR_ij = 1.0 / r_ij;
      pre_thrInvR_ij[neighbor_j] = invR_ij;

      invRMinusCutoffA = 1.0 / (r_ij - cutoffA);
      sigmaInvRMinusCutoffA = sigma * invRMinusCutoffA;
      gammInvRMinusCutoffA = gamm * invRMinusCutoffA;

      interpolDeltaX = r_ij - GRIDSTART;
      interpolTMP = (interpolDeltaX * GRIDDENSITY);
      interpolIDX = (int) interpolTMP;

      interpolY1 = exp3B[interpolIDX];
      interpolY2 = exp3B[interpolIDX + 1];
      exp3B_ij = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

      exp3BDerived_ij = -exp3B_ij * gammInvRMinusCutoffA * invRMinusCutoffA;

      pre_thrExp3B_ij[neighbor_j] = exp3B_ij;
      pre_thrExp3BDerived_ij[neighbor_j] = exp3BDerived_ij;

      interpolY1 = exp2B[interpolIDX];
      interpolY2 = exp2B[interpolIDX + 1];
      exp2B_ij = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

      exp2BDerived_ij = -exp2B_ij * sigmaInvRMinusCutoffA * invRMinusCutoffA;

      pre_thrExp2B_ij[neighbor_j] = exp2B_ij;
      pre_thrExp2BDerived_ij[neighbor_j] = exp2BDerived_ij;

      interpolY1 = pow2B[interpolIDX];
      interpolY2 = pow2B[interpolIDX + 1];
      pow2B_ij = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

      pre_thrPow2B_ij[neighbor_j] = pow2B_ij;

      // zeta and its derivative

      if (r_ij < cutoffC)
        zeta_i += 1.0;
      else {
        interpolY1 = cutoffFunction[interpolIDX];
        interpolY2 = cutoffFunction[interpolIDX + 1];
        cutoffFunction_ij = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

        zeta_i += cutoffFunction_ij;

        interpolY1 = cutoffFunctionDerived[interpolIDX];
        interpolY2 = cutoffFunctionDerived[interpolIDX + 1];
        zeta_iDerived = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

        zeta_iDerivedInvR_ij = zeta_iDerived * invR_ij;

        preForceCoord_counter = numForceCoordPairs * 5;
        pre_thrForceCoord[preForceCoord_counter + 0] = zeta_iDerivedInvR_ij;
        pre_thrForceCoord[preForceCoord_counter + 1] = dr_ij[0];
        pre_thrForceCoord[preForceCoord_counter + 2] = dr_ij[1];
        pre_thrForceCoord[preForceCoord_counter + 3] = dr_ij[2];
        pre_thrForceCoord[preForceCoord_counter + 4] = j;
        numForceCoordPairs++;
      }
    }

    // quantities depending on zeta_i

    interpolDeltaX = zeta_i;
    interpolTMP = (interpolDeltaX * GRIDDENSITY);
    interpolIDX = (int) interpolTMP;

    interpolY1 = expMinusBetaZeta_iZeta_iGrid[interpolIDX];
    interpolY2 = expMinusBetaZeta_iZeta_iGrid[interpolIDX + 1];
    expMinusBetaZeta_iZeta_i = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

    interpolY1 = qFunctionGrid[interpolIDX];
    interpolY2 = qFunctionGrid[interpolIDX + 1];
    qFunction = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

    interpolY1 = tauFunctionGrid[interpolIDX];
    interpolY2 = tauFunctionGrid[interpolIDX + 1];
    tauFunction = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

    interpolY1 = tauFunctionDerivedGrid[interpolIDX];
    interpolY2 = tauFunctionDerivedGrid[interpolIDX + 1];
    tauFunctionDerived = interpolY1 + (interpolY2 - interpolY1) * (interpolTMP - interpolIDX);

    forceModCoord_factor = 2.0 * beta * zeta_i * expMinusBetaZeta_iZeta_i;

    forceModCoord = 0.0;

    // two-body interactions, skip half of them

    for (int neighbor_j = 0; neighbor_j < jnum; neighbor_j++) {
      double dr_ij[3], r_ij, f_ij[3];

      j = jlist[neighbor_j];
      j &= NEIGHMASK;

      dr_ij[0] = x[j].x - xtmp;
      dr_ij[1] = x[j].y - ytmp;
      dr_ij[2] = x[j].z - ztmp;
      r_ij = dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1] + dr_ij[2] * dr_ij[2];

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (r_ij > params[ijparam].cutsq) continue;

      r_ij = sqrt(r_ij);

      invR_ij = pre_thrInvR_ij[neighbor_j];
      pow2B_ij = pre_thrPow2B_ij[neighbor_j];

      potential2B_factor = pow2B_ij - expMinusBetaZeta_iZeta_i;

      exp2B_ij = pre_thrExp2B_ij[neighbor_j];

      pow2BDerived_ij = -rho * invR_ij * pow2B_ij;

      forceModCoord += (forceModCoord_factor * exp2B_ij);

      exp2BDerived_ij = pre_thrExp2BDerived_ij[neighbor_j];
      forceMod2B = exp2BDerived_ij * potential2B_factor + exp2B_ij * pow2BDerived_ij;

      directorCos_ij_x = invR_ij * dr_ij[0];
      directorCos_ij_y = invR_ij * dr_ij[1];
      directorCos_ij_z = invR_ij * dr_ij[2];

      exp3B_ij = pre_thrExp3B_ij[neighbor_j];
      exp3BDerived_ij = pre_thrExp3BDerived_ij[neighbor_j];

      f_ij[0] = forceMod2B * directorCos_ij_x;
      f_ij[1] = forceMod2B * directorCos_ij_y;
      f_ij[2] = forceMod2B * directorCos_ij_z;

      f[i].x += f_ij[0];
      f[i].y += f_ij[1];
      f[i].z += f_ij[2];

      f[j].x -= f_ij[0];
      f[j].y -= f_ij[1];
      f[j].z -= f_ij[2];

      // potential energy

      evdwl = (exp2B_ij * potential2B_factor);

      if (EVFLAG)
        ev_tally_thr(this, i, j, nlocal, /* newton_pair */ 1, evdwl, 0.0, -forceMod2B * invR_ij,
                     dr_ij[0], dr_ij[1], dr_ij[2], thr);

      // three-body Forces

      for (int neighbor_k = neighbor_j + 1; neighbor_k < jnum; neighbor_k++) {
        double dr_ik[3], r_ik, f_ik[3];

        k = jlist[neighbor_k];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];

        dr_ik[0] = x[k].x - xtmp;
        dr_ik[1] = x[k].y - ytmp;
        dr_ik[2] = x[k].z - ztmp;
        r_ik = dr_ik[0] * dr_ik[0] + dr_ik[1] * dr_ik[1] + dr_ik[2] * dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        r_ik = sqrt(r_ik);

        invR_ik = pre_thrInvR_ij[neighbor_k];

        directorCos_ik_x = invR_ik * dr_ik[0];
        directorCos_ik_y = invR_ik * dr_ik[1];
        directorCos_ik_z = invR_ik * dr_ik[2];

        cosTeta = directorCos_ij_x * directorCos_ik_x + directorCos_ij_y * directorCos_ik_y +
            directorCos_ij_z * directorCos_ik_z;

        cosTetaDiff = cosTeta + tauFunction;
        cosTetaDiffCosTetaDiff = cosTetaDiff * cosTetaDiff;
        qFunctionCosTetaDiffCosTetaDiff = cosTetaDiffCosTetaDiff * qFunction;
        expMinusQFunctionCosTetaDiffCosTetaDiff = exp(-qFunctionCosTetaDiffCosTetaDiff);

        potentia3B_factor = lambda *
            ((1.0 - expMinusQFunctionCosTetaDiffCosTetaDiff) +
             eta * qFunctionCosTetaDiffCosTetaDiff);

        exp3B_ik = pre_thrExp3B_ij[neighbor_k];
        exp3BDerived_ik = pre_thrExp3BDerived_ij[neighbor_k];

        forceMod3B_factor1_ij = -exp3BDerived_ij * exp3B_ik * potentia3B_factor;
        forceMod3B_factor2 = 2.0 * lambda * exp3B_ij * exp3B_ik * qFunction * cosTetaDiff *
            (eta + expMinusQFunctionCosTetaDiffCosTetaDiff);
        forceMod3B_factor2_ij = forceMod3B_factor2 * invR_ij;

        f_ij[0] = forceMod3B_factor1_ij * directorCos_ij_x +
            forceMod3B_factor2_ij * (cosTeta * directorCos_ij_x - directorCos_ik_x);
        f_ij[1] = forceMod3B_factor1_ij * directorCos_ij_y +
            forceMod3B_factor2_ij * (cosTeta * directorCos_ij_y - directorCos_ik_y);
        f_ij[2] = forceMod3B_factor1_ij * directorCos_ij_z +
            forceMod3B_factor2_ij * (cosTeta * directorCos_ij_z - directorCos_ik_z);

        forceMod3B_factor1_ik = -exp3BDerived_ik * exp3B_ij * potentia3B_factor;
        forceMod3B_factor2_ik = forceMod3B_factor2 * invR_ik;

        f_ik[0] = forceMod3B_factor1_ik * directorCos_ik_x +
            forceMod3B_factor2_ik * (cosTeta * directorCos_ik_x - directorCos_ij_x);
        f_ik[1] = forceMod3B_factor1_ik * directorCos_ik_y +
            forceMod3B_factor2_ik * (cosTeta * directorCos_ik_y - directorCos_ij_y);
        f_ik[2] = forceMod3B_factor1_ik * directorCos_ik_z +
            forceMod3B_factor2_ik * (cosTeta * directorCos_ik_z - directorCos_ij_z);

        forceModCoord += (forceMod3B_factor2 * (tauFunctionDerived - 0.5 * mu * cosTetaDiff));

        f[j].x += f_ij[0];
        f[j].y += f_ij[1];
        f[j].z += f_ij[2];

        f[k].x += f_ik[0];
        f[k].y += f_ik[1];
        f[k].z += f_ik[2];

        f[i].x -= f_ij[0] + f_ik[0];
        f[i].y -= f_ij[1] + f_ik[1];
        f[i].z -= f_ij[2] + f_ik[2];

        // potential energy

        evdwl = (exp3B_ij * exp3B_ik * potentia3B_factor);

        if (EVFLAG) ev_tally3_thr(this, i, j, k, evdwl, 0.0, f_ij, f_ik, dr_ij, dr_ik, thr);
      }
    }

    // forces due to environment coordination f(Z)

    for (int idx = 0; idx < numForceCoordPairs; idx++) {
      double dr_ij[3], f_ij[3];

      preForceCoord_counter = idx * 5;
      zeta_iDerivedInvR_ij = pre_thrForceCoord[preForceCoord_counter + 0];
      dr_ij[0] = pre_thrForceCoord[preForceCoord_counter + 1];
      dr_ij[1] = pre_thrForceCoord[preForceCoord_counter + 2];
      dr_ij[2] = pre_thrForceCoord[preForceCoord_counter + 3];
      j = static_cast<int>(pre_thrForceCoord[preForceCoord_counter + 4]);

      forceModCoord_ij = forceModCoord * zeta_iDerivedInvR_ij;

      f_ij[0] = forceModCoord_ij * dr_ij[0];
      f_ij[1] = forceModCoord_ij * dr_ij[1];
      f_ij[2] = forceModCoord_ij * dr_ij[2];

      f[i].x -= f_ij[0];
      f[i].y -= f_ij[1];
      f[i].z -= f_ij[2];

      f[j].x += f_ij[0];
      f[j].y += f_ij[1];
      f[j].z += f_ij[2];

      if (EVFLAG)
        ev_tally_thr(this, i, j, nlocal, /* newton_pair */ 1, 0.0, 0.0, -forceModCoord_ij, dr_ij[0],
                     dr_ij[1], dr_ij[2], thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairEDIPOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairEDIP::memory_usage();

  return bytes;
}

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "pair_lj_cut_coul_cut_dielectric_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 1e-6

/* ---------------------------------------------------------------------- */

PairLJCutCoulCutDielectricOMP::PairLJCutCoulCutDielectricOMP(LAMMPS *lmp) :
    PairLJCutCoulCutDielectric(lmp), ThrOMP(lmp, THR_PAIR)
{
}

/* ---------------------------------------------------------------------- */

PairLJCutCoulCutDielectricOMP::~PairLJCutCoulCutDielectricOMP() {}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulCutDielectricOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (atom->nmax > nmax) {
    memory->destroy(efield);
    memory->destroy(epot);
    nmax = atom->nmax;
    memory->create(efield, nmax, 3, "pair:efield");
    memory->create(epot, nmax, "pair:epot");
  }

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
        if (force->newton_pair)
          eval<1, 1, 1>(ifrom, ito, thr);
        else
          eval<1, 1, 0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair)
          eval<1, 0, 1>(ifrom, ito, thr);
        else
          eval<1, 0, 0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair)
        eval<0, 0, 1>(ifrom, ito, thr);
      else
        eval<0, 0, 0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJCutCoulCutDielectricOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double qtmp, etmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul;
  double fpair_i, fpair_j;
  double rsq, r2inv, r6inv, forcecoul, forcelj, factor_coul, factor_lj;
  double efield_i, epot_i;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = ecoul = 0.0;

  const dbl3_t *_noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const double *_noalias const q = atom->q;
  const double *_noalias const eps = atom->epsilon;
  const dbl3_t *_noalias const norm = (dbl3_t *) atom->mu[0];
  const double *_noalias const curvature = atom->curvature;
  const double *_noalias const area = atom->area;
  const int *_noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double *_noalias const special_coul = force->special_coul;
  const double *_noalias const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  double fxtmp, fytmp, fztmp, extmp, eytmp, eztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    etmp = eps[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;
    extmp = eytmp = eztmp = 0.0;

    // self term Eq. (55) for I_{ii} and Eq. (52) and in Barros et al
    double curvature_threshold = sqrt(area[i]);
    if (curvature[i] < curvature_threshold) {
      double sf = curvature[i] / (4.0 * MY_PIS * curvature_threshold) * area[i] * q[i];
      efield[i][0] = sf * norm[i].x;
      efield[i][1] = sf * norm[i].y;
      efield[i][2] = sf * norm[i].z;
    } else {
      efield[i][0] = efield[i][1] = efield[i][2] = 0;
    }

    epot[i] = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0 / rsq;

        if (rsq < cut_coulsq[itype][jtype] && rsq > EPSILON) {
          efield_i = q[j] * sqrt(r2inv);
          forcecoul = qqrd2e * qtmp * efield_i;
          epot_i = efield_i;
        } else
          epot_i = efield_i = forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
        } else
          forcelj = 0.0;

        fpair_i = (factor_coul * etmp * forcecoul + factor_lj * forcelj) * r2inv;

        fxtmp += delx * fpair_i;
        fytmp += dely * fpair_i;
        fztmp += delz * fpair_i;

        efield_i *= (factor_coul * etmp * r2inv);
        extmp += delx * efield_i;
        eytmp += dely * efield_i;
        eztmp += delz * efield_i;
        epot[i] += epot_i;

        if (NEWTON_PAIR || j >= nlocal) {
          fpair_j = (factor_coul * eps[j] * forcecoul + factor_lj * forcelj) * r2inv;
          f[j].x -= delx * fpair_j;
          f[j].y -= dely * fpair_j;
          f[j].z -= delz * fpair_j;
        }

        if (EFLAG) {
          if (rsq < cut_coulsq[itype][jtype]) {
            ecoul = factor_coul * qqrd2e * qtmp * q[j] * (etmp + eps[j]) * sqrt(r2inv);
          } else
            ecoul = 0.0;
          ecoul *= 0.5;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
            evdwl *= factor_lj;
          } else
            evdwl = 0.0;
        }

        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, ecoul, fpair_i, delx, dely, delz,
                       thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
    efield[i][0] += extmp;
    efield[i][1] += eytmp;
    efield[i][2] += eztmp;
  }
}

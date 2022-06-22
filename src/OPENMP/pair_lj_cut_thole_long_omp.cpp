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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "pair_lj_cut_thole_long_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_drude.h"
#include "force.h"
#include "math_const.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F   1.12837917
#define EWALD_P   9.95473818e-1
#define B0       -0.1335096380159268
#define B1       -2.57839507e-1
#define B2       -1.37203639e-1
#define B3       -8.88822059e-3
#define B4       -5.80844129e-3
#define B5        1.14652755e-1

#define EPSILON 1.0e-20
#define EPS_EWALD 1.0e-6
#define EPS_EWALD_SQR 1.0e-12

/* ---------------------------------------------------------------------- */

PairLJCutTholeLongOMP::PairLJCutTholeLongOMP(LAMMPS *lmp) :
    PairLJCutTholeLong(lmp), ThrOMP(lmp, THR_PAIR)
{
    suffix_flag |= Suffix::OMP;
    respa_enable = 0;
    cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairLJCutTholeLongOMP::compute(int eflag, int vflag)
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

    if (evflag) {
      if (eflag) {
        if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }
    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJCutTholeLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * const q = atom->q;
  const int * _noalias const type = atom->type;
  const double * _noalias const special_lj = force->special_lj;
  const double * _noalias const special_coul = force->special_coul;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;
  const int * _noalias const drudetype = fix_drude->drudetype;
  const tagint * _noalias const drudeid = fix_drude->drudeid;

  double xtmp,ytmp,ztmp,delx,dely,delz,fxtmp,fytmp,fztmp;

  const int nlocal = atom->nlocal;

  int j,jj,jnum,jtype,itable;
  double ecoul,fpair,evdwl;
  double r,rsq,r2inv,forcecoul,factor_coul,forcelj,factor_lj,r6inv;
  double fraction,table;
  double grij,expm2,prefactor,t,erfc,u;
  double factor_f,factor_e;
  int di,dj;
  double qj,dqi,dqj,dcoul,asr,exp_asr;
  int di_closest;
  const double qqrd2e = force->qqrd2e;

  evdwl = ecoul = 0.0;

  // loop over neighbors of my atoms
  for (int ii = iifrom; ii < iito; ii++) {
    const int i = ilist[ii];
    const double qi = q[i];
    const int itype = type[i];
    const int * _noalias const jlist = firstneigh[i];
    const double * _noalias const cutsqi = cutsq[itype];
    const double * _noalias const cut_ljsqi = cut_ljsq[itype];
    const double * _noalias const offseti = offset[itype];
    const double * _noalias const lj1i = lj1[itype];
    const double * _noalias const lj2i = lj2[itype];
    const double * _noalias const lj3i = lj3[itype];
    const double * _noalias const lj4i = lj4[itype];

    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.;

    if (drudetype[type[i]] != NOPOL_TYPE) {
      di = atom->map(drudeid[i]);
      if (di < 0) error->all(FLERR, "Drude partner not found");
      di_closest = domain->closest_image(i, di);
      if (drudetype[type[i]] == CORE_TYPE)
        dqi = -q[di];
      else
        dqi = qi;
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsqi[jtype]) {
        rsq += EPSILON; // Add Epsilon for case: r = 0; DC-DP 1-1 interaction must be removed by special bond;
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          qj = q[j];
          r = sqrt(rsq);

          if (!ncoultablebits || rsq <= tabinnersq) {
            grij = g_ewald * (r + EPS_EWALD);
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            u = 1. - t;
            erfc = t * (1.+u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
            prefactor = qqrd2e * qi*qj/(r + EPS_EWALD);
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
            r2inv = 1.0/(rsq + EPS_EWALD_SQR);
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qi*qj * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qi*qj * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }

          if (drudetype[type[i]] != NOPOL_TYPE &&
              drudetype[type[j]] != NOPOL_TYPE) {
            if (j != di_closest) {
              if (drudetype[type[j]] == CORE_TYPE) {
                dj = atom->map(drudeid[j]);
                dqj = -q[dj];
              } else dqj = qj;
              asr = ascreen[type[i]][type[j]] * r;
              exp_asr = exp(-asr);
              dcoul = qqrd2e * dqi * dqj / r;
              factor_f = 0.5*(2. + (exp_asr * (-2. - asr * (2. + asr))))
                  - factor_coul;
              if (EFLAG) factor_e = 0.5*(2. - (exp_asr * (2. + asr)))
                             - factor_coul;
              forcecoul += factor_f * dcoul;
            }
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsqi[jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1i[jtype]*r6inv - lj2i[jtype]);
        } else forcelj = 0.0;

        fpair = (forcecoul + factor_lj*forcelj) * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG) {
          if (rsq < cut_coulsq) {
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor*erfc;
            else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qi*qj * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
            if (drudetype[type[i]] != NOPOL_TYPE &&
                drudetype[type[j]] != NOPOL_TYPE && j != di_closest) {
              ecoul += factor_e * dcoul;
            }
          } else ecoul = 0.0;

          if (rsq < cut_ljsqi[jtype]) {
            evdwl = r6inv*(lj3i[jtype]*r6inv-lj4i[jtype]) -
              offseti[jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (EVFLAG) ev_tally_thr(this, i,j,nlocal,NEWTON_PAIR,
                                 evdwl,ecoul,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}


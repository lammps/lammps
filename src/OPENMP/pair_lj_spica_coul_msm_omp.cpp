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
   This style is a simplified re-implementation of the CG/CMM pair style
------------------------------------------------------------------------- */

#include "pair_lj_spica_coul_msm_omp.h"
#include "lj_spica_common.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace LJSPICAParms;
/* ---------------------------------------------------------------------- */

PairLJSPICACoulMSMOMP::PairLJSPICACoulMSMOMP(LAMMPS *lmp) :
  PairLJSPICACoulMSM(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJSPICACoulMSMOMP::compute(int eflag, int vflag)
{
  if (force->kspace->scalar_pressure_flag)
    error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' "
      "with OMP MSM Pair styles");

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
        if (force->newton_pair) eval_msm_thr<1,1,1>(ifrom, ito, thr);
        else eval_msm_thr<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval_msm_thr<1,0,1>(ifrom, ito, thr);
        else eval_msm_thr<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval_msm_thr<0,0,1>(ifrom, ito, thr);
      else eval_msm_thr<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJSPICACoulMSMOMP::eval_msm_thr(int iifrom, int iito, ThrData * const thr)
{

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;

  const int * const ilist = list->ilist;
  const int * const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;
  const int nlocal = atom->nlocal;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {

    const int i = ilist[ii];
    const int itype = type[i];
    const double qtmp = q[i];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    double fxtmp,fytmp,fztmp;
    fxtmp=fytmp=fztmp=0.0;

    const int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      double forcecoul, forcelj, evdwl, ecoul, fgamma, egamma;
      forcecoul = forcelj = evdwl = ecoul = 0.0;

      const int sbindex = sbmask(jlist[jj]);
      const int j = jlist[jj] & NEIGHMASK;

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;
      const int jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        const double r2inv = 1.0/rsq;
        const int ljt = lj_type[itype][jtype];

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            const double r = sqrt(rsq);
            const double prefactor = qqrd2e * qtmp*q[j]/r;
            fgamma = 1.0 + (rsq/cut_coulsq)*force->kspace->dgamma(r/cut_coul);
            forcecoul = prefactor * fgamma;
            if (EFLAG) {
              egamma = 1.0 - (r/cut_coul)*force->kspace->gamma(r/cut_coul);
              ecoul = prefactor*egamma;
            }
            if (sbindex) {
              const double adjust = (1.0-special_coul[sbindex])*prefactor;
              forcecoul -= adjust;
              if (EFLAG) ecoul -= adjust;
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
            const double fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            const double table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (EFLAG) ecoul = qtmp*q[j] * (etable[itable] + fraction*detable[itable]);
            if (sbindex) {
              const double table2 = ctable[itable] + fraction*dctable[itable];
              const double prefactor = qtmp*q[j] * table2;
              const double adjust = (1.0-special_coul[sbindex])*prefactor;
              forcecoul -= adjust;
              if (EFLAG) ecoul -= adjust;
            }
          }
        }

        if (rsq < cut_ljsq[itype][jtype]) {

          if (ljt == LJ12_4) {
            const double r4inv=r2inv*r2inv;
            forcelj = r4inv*(lj1[itype][jtype]*r4inv*r4inv
                             - lj2[itype][jtype]);

            if (EFLAG)
              evdwl = r4inv*(lj3[itype][jtype]*r4inv*r4inv
                             - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ9_6) {
            const double r3inv = r2inv*sqrt(r2inv);
            const double r6inv = r3inv*r3inv;
            forcelj = r6inv*(lj1[itype][jtype]*r3inv
                             - lj2[itype][jtype]);
            if (EFLAG)
              evdwl = r6inv*(lj3[itype][jtype]*r3inv
                             - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ12_6) {
            const double r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv*(lj1[itype][jtype]*r6inv
                             - lj2[itype][jtype]);
            if (EFLAG)
              evdwl = r6inv*(lj3[itype][jtype]*r6inv
                             - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ12_5) {
            const double r5inv = r2inv*r2inv*sqrt(r2inv);
            const double r7inv = r5inv*r2inv;
            forcelj = r5inv*(lj1[itype][jtype]*r7inv
                             - lj2[itype][jtype]);
            if (EFLAG)
              evdwl = r5inv*(lj3[itype][jtype]*r7inv
                             - lj4[itype][jtype]) - offset[itype][jtype];
          }

          if (sbindex) {
            const double factor_lj = special_lj[sbindex];
            forcelj *= factor_lj;
            if (EFLAG) evdwl *= factor_lj;
          }

        }

        const double fpair = (forcecoul + forcelj) * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 evdwl,ecoul,fpair,delx,dely,delz,thr);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJSPICACoulMSMOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJSPICACoulMSM::memory_usage();

  return bytes;
}

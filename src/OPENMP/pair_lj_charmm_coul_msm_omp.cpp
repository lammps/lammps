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

#include "pair_lj_charmm_coul_msm_omp.h"

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

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulMSMOMP::PairLJCharmmCoulMSMOMP(LAMMPS *lmp) :
  PairLJCharmmCoulMSM(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulMSMOMP::compute(int eflag, int vflag)
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

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJCharmmCoulMSMOMP::eval(int iifrom, int iito, ThrData * const thr)
{

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const double * _noalias const special_coul = force->special_coul;
  const double * _noalias const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  const double inv_denom_lj = 1.0/denom_lj;

  const int * const ilist = list->ilist;
  const int * const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;
  const int nlocal = atom->nlocal;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {

    const int i = ilist[ii];
    const int itype = type[i];
    const double qtmp = q[i];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    double fxtmp,fytmp,fztmp;
    fxtmp=fytmp=fztmp=0.0;

    const int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      double forcecoul, forcelj, evdwl, ecoul;
      forcecoul = forcelj = evdwl = ecoul = 0.0;

      const int sbindex = sbmask(jlist[jj]);
      const int j = jlist[jj] & NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx*delx + dely*dely + delz*delz;
      const int jtype = type[j];

      if (rsq < cut_bothsq) {
        const double r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {

            const double r = sqrt(rsq);
            const double prefactor = qqrd2e * qtmp*q[j]/r;
            const double egamma = 1.0 - (r/cut_coul)*force->kspace->gamma(r/cut_coul);
            const double fgamma = 1.0 + (rsq/cut_coulsq)*force->kspace->dgamma(r/cut_coul);
            forcecoul = prefactor * fgamma;

            if (EFLAG) ecoul = prefactor*egamma;
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

        if (rsq < cut_ljsq) {
          const double r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          if (EFLAG) evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);

          if (rsq > cut_lj_innersq) {
            const double drsq = cut_ljsq - rsq;
            const double cut2 = (rsq - cut_lj_innersq) * drsq;
            const double switch1 = drsq * (drsq*drsq + 3.0*cut2) * inv_denom_lj;
            const double switch2 = 12.0*rsq * cut2 * inv_denom_lj;
            if (EFLAG) {
              forcelj = forcelj*switch1 + evdwl*switch2;
              evdwl *= switch1;
            } else {
              const double philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
              forcelj =  forcelj*switch1 + philj*switch2;
            }
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
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 evdwl,ecoul,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCharmmCoulMSMOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJCharmmCoulMSM::memory_usage();

  return bytes;
}

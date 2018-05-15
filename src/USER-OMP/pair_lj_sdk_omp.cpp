/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
   This style is a simplified re-implementation of the CG/CMM pair style
------------------------------------------------------------------------- */

#include <cmath>
#include "pair_lj_sdk_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "lj_sdk_common.h"

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace LJSDKParms;

/* ---------------------------------------------------------------------- */

PairLJSDKOMP::PairLJSDKOMP(LAMMPS *lmp) :
  PairLJSDK(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJSDKOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (evflag) {
      if (eflag) {
        if (force->newton_pair) eval_thr<1,1,1>(ifrom, ito, thr);
        else eval_thr<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval_thr<1,0,1>(ifrom, ito, thr);
        else eval_thr<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval_thr<0,0,1>(ifrom, ito, thr);
      else eval_thr<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJSDKOMP::eval_thr(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,forcelj,factor_lj;

  evdwl = 0.0;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  double fxtmp,fytmp,fztmp;

  const int * const ilist = list->ilist;
  const int * const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    fxtmp=fytmp=fztmp=0.0;

    const int itype = type[i];
    const int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        const int ljt = lj_type[itype][jtype];

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
        } else continue;

        fpair = factor_lj*forcelj*r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG) evdwl *= factor_lj;
        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 evdwl,0.0,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJSDKOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJSDK::memory_usage();

  return bytes;
}

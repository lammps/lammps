/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include <cmath>
#include "pair_nm_cut_coul_cut_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairNMCutCoulCutOMP::PairNMCutCoulCutOMP(LAMMPS *lmp) :
  PairNMCutCoulCut(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairNMCutCoulCutOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

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
void PairNMCutCoulCutOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int j,ii,jj,jnum,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r,rsq,r2inv,rminv,rninv,forcecoul,forcenm,factor_coul,factor_lj;
  int *ilist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const double * _noalias const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    const int i = ilist[ii];
    const int itype = type[i];
    const int    * _noalias const jlist = firstneigh[i];
    const double * _noalias const cutsqi = cutsq[itype];
    const double * _noalias const cut_coulsqi = cut_coulsq[itype];
    const double * _noalias const cut_ljsqi = cut_ljsq[itype];
    const double * _noalias const offseti = offset[itype];
    const double * _noalias const mmi = mm[itype];
    const double * _noalias const nni = nn[itype];
    const double * _noalias const nmi = nm[itype];
    const double * _noalias const e0nmi = e0nm[itype];
    const double * _noalias const r0mi = r0m[itype];
    const double * _noalias const r0ni = r0n[itype];

    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

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
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsqi[jtype]) {
          const double rinv = sqrt(r2inv);
          forcecoul = qqrd2e * qtmp*q[j]*rinv;
          forcecoul *= factor_coul;
          if (EFLAG) ecoul = factor_coul * qqrd2e * qtmp*q[j]*rinv;
        } else {
          forcecoul = 0.0;
          if (EFLAG) ecoul = 0.0;
        }

        if (rsq < cut_ljsqi[jtype]) {
          r = sqrt(rsq);
          rminv = pow(r2inv,mmi[jtype]/2.0);
          rninv = pow(r2inv,nni[jtype]/2.0);
          forcenm = e0nmi[jtype]*nmi[jtype] *
            (r0ni[jtype]/pow(r,nni[jtype]) -
             r0mi[jtype]/pow(r,mmi[jtype]));
          forcenm *= factor_lj;
          if (EFLAG)
            evdwl = (e0nmi[jtype]*(mmi[jtype] *
                                   r0ni[jtype]*rninv -
                                   nni[jtype] *
                                   r0mi[jtype]*rminv) -
                     offseti[jtype]) * factor_lj;
        } else {
          forcenm = 0.0;
          if (EFLAG) evdwl = 0.0;
        }

        fpair = (forcecoul + forcenm) * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
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

/* ---------------------------------------------------------------------- */

double PairNMCutCoulCutOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairNMCutCoulCut::memory_usage();

  return bytes;
}

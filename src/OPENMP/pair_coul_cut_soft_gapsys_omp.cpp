// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_coul_cut_soft_gapsys_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulCutSoftGapsysOMP::PairCoulCutSoftGapsysOMP(LAMMPS *lmp) :
  PairCoulCutSoftGapsys(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairCoulCutSoftGapsysOMP::compute(int eflag, int vflag)
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

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairCoulCutSoftGapsysOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,forcecoul,factor_coul;
  double cut_inner;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double qqrd2e = force->qqrd2e;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      cut_inner = (1.0 + sigmaq * fabs(qtmp * q[j])) * alphaq * pow(lambda[itype][jtype], 1.0 / 6.0);

      if (rsq < cut_inner * cut_inner) {

        forcecoul = factor_coul * qqrd2e * qtmp * q[j];
        fpair = (- 2.0 / pow(cut_inner, 3) + 3.0 / (pow(cut_inner, 2) * sqrt(rsq))) * forcecoul;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (EFLAG)
          ecoul = factor_coul * qqrd2e * qtmp * q[j] *
          (rsq / pow(cut_inner, 3) - 3 * sqrt(rsq) / pow(cut_inner, 2) + 3 / cut_inner);

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 0.0,ecoul,fpair,delx,dely,delz,thr);

      } else if (rsq < cutsq[itype][jtype]) {

        double r2inv = 1.0 / rsq;
        double rinv = sqrt(r2inv);
        forcecoul = qqrd2e * qtmp * q[j] * rinv;
        fpair = factor_coul * forcecoul * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (EFLAG)
          ecoul = factor_coul * qqrd2e * qtmp * q[j] * rinv;

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 0.0,ecoul,fpair,delx,dely,delz,thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairCoulCutSoftGapsysOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairCoulCutSoftGapsys::memory_usage();

  return bytes;
}

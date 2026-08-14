// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_nm_cut_split_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_special.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;
using MathSpecial::powint;

/* ---------------------------------------------------------------------- */

PairNMCutSplitOMP::PairNMCutSplitOMP(LAMMPS *lmp) :
  PairNMCutSplit(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairNMCutSplitOMP::compute(int eflag, int vflag)
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
void PairNMCutSplitOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,factor_lj;
  double r,forcenm,rminv,rninv;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_lj = force->special_lj;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);

        // r < r0 --> use generalized LJ
        if (rsq < r0[itype][jtype]*r0[itype][jtype]) {
          forcenm = e0nm[itype][jtype]*nm[itype][jtype]*
            (r0n[itype][jtype]/pow(r,nn[itype][jtype])
             -r0m[itype][jtype]/pow(r,mm[itype][jtype]));
        }
        // r > r0 --> use standard LJ (m = 6 n = 12)
        else forcenm =(e0[itype][jtype]/6.0)*72.0*(4.0/powint(r,12)-2.0/powint(r,6));

        fpair = factor_lj*forcenm*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (EFLAG) {
          // r < r0 --> use generalized LJ
          if (rsq < r0[itype][jtype]*r0[itype][jtype]) {
            rminv = pow(r2inv,mm[itype][jtype]/2.0);
            rninv = pow(r2inv,nn[itype][jtype]/2.0);

            evdwl = e0nm[itype][jtype]*(mm[itype][jtype]*r0n[itype][jtype]*rninv -
                                        nn[itype][jtype]*r0m[itype][jtype]*rminv) -
              offset[itype][jtype];
          }
          // r > r0 --> use standard LJ (m = 6 n = 12)
          else evdwl = (e0[itype][jtype]/6.0)*(24.0*powint(r2inv,6) - 24.0*pow(r2inv,3));
          evdwl *= factor_lj;
        }
        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,evdwl,0.0,fpair,delx,dely,delz,thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairNMCutSplitOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairNMCutSplit::memory_usage();

  return bytes;
}

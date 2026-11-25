// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Chun-Yeol You (DGIST)
------------------------------------------------------------------------- */

#include "pair_rheo_solid_omp.h"

#include "atom.h"
#include "comm.h"
#include "fix_rheo.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace RHEO_NS;

/* ---------------------------------------------------------------------- */

PairRHEOSolidOMP::PairRHEOSolidOMP(LAMMPS *lmp) :
  PairRHEOSolid(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairRHEOSolidOMP::compute(int eflag, int vflag)
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
        if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else                    eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else                    eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else                    eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairRHEOSolidOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  double * const * const v    = atom->v;
  const int * _noalias const type   = atom->type;
  const int * _noalias const status = atom->rheo_status;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;

  const int * _noalias const ilist    = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  int ** const firstneigh             = list->firstneigh;

  for (int ii = iifrom; ii < iito; ii++) {
    const int i = ilist[ii];
    if (!(status[i] & STATUS_SOLID)) continue;

    const double xtmp  = x[i].x;
    const double ytmp  = x[i].y;
    const double ztmp  = x[i].z;
    const double vxtmp = v[i][0];
    const double vytmp = v[i][1];
    const double vztmp = v[i][2];
    const int itype    = type[i];
    const int * _noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    double fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_lj = special_lj[sbmask(j)];
      if (factor_lj == 0) continue;
      j &= NEIGHMASK;

      if (!(status[j] & STATUS_SOLID)) continue;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq  = delx*delx + dely*dely + delz*delz;
      const int jtype   = type[j];

      if (rsq < cutsq[itype][jtype]) {
        const double r   = sqrt(rsq);
        const double rinv = 1.0/r;

        double fpair = k[itype][jtype] * (cut[itype][jtype] - r);

        double smooth = rsq / cutsq[itype][jtype];
        smooth *= smooth;
        smooth *= smooth;
        smooth = 1.0 - smooth;

        const double delvx = vxtmp - v[j][0];
        const double delvy = vytmp - v[j][1];
        const double delvz = vztmp - v[j][2];
        const double dot = delx*delvx + dely*delvy + delz*delvz;
        fpair -= gamma[itype][jtype] * dot * smooth * rinv;

        fpair *= factor_lj * rinv;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        const double evdwl = 0.0;
        if (EVFLAG) ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR,
                                  evdwl, 0.0, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairRHEOSolidOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairRHEOSolid::memory_usage();
  return bytes;
}

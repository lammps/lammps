// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "omp_compat.h"
#include "pair_atm_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairATMOMP::PairATMOMP(LAMMPS *lmp) :
  PairATM(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairATMOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

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
      if (eflag) eval<1,1>(ifrom, ito, thr);
      else eval<1,0>(ifrom, ito, thr);
    } else {
      eval<0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG>
void PairATMOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i, j, k, ii, jj, kk, jnum, jnumm1;
  double xi, yi, zi, evdwl;
  double rij2, rik2, rjk2;
  double rij[3], rik[3], rjk[3], fj[3], fk[3];
  double nu_local;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;

  const double cutoff_squared = cut_global * cut_global;
  const double triple = cut_triple * cut_triple * cut_triple;
  const double cutoff_triple_sixth = triple * triple;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // triple loop over local atoms and neighbors twice
  // must compute each IJK triplet interaction exactly once
  // newton_pair is always on for ATM; uses full neighbor list

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    xi = x[i].x;
    yi = x[i].y;
    zi = x[i].z;

    jlist = firstneigh[i];
    jnum = numneigh[i];
    jnumm1 = jnum - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      rij[0] = x[j].x - xi;
      if (rij[0] < 0.0) continue;
      rij[1] = x[j].y - yi;
      if (rij[0] == 0.0 && rij[1] < 0.0) continue;
      rij[2] = x[j].z - zi;
      if (rij[0] == 0.0 && rij[1] == 0.0 && rij[2] < 0.0) continue;
      rij2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
      if (rij2 > cutoff_squared) continue;

      for (kk = jj+1; kk < jnum; kk++) {
        k = jlist[kk];
        k &= NEIGHMASK;

        rik[0] = x[k].x - xi;
        if (rik[0] < 0.0) continue;
        rik[1] = x[k].y - yi;
        if (rik[0] == 0.0 && rik[1] < 0.0) continue;
        rik[2] = x[k].z - zi;
        if (rik[0] == 0.0 && rik[1] == 0.0 && rik[2] < 0.0) continue;
        rik2 = rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2];
        if (rik2 > cutoff_squared) continue;

        rjk[0] = x[k].x - x[j].x;
        rjk[1] = x[k].y - x[j].y;
        rjk[2] = x[k].z - x[j].z;
        rjk2 = rjk[0]*rjk[0] + rjk[1]*rjk[1] + rjk[2]*rjk[2];
        if (rjk2 > cutoff_squared) continue;

        double r6 = rij2 * rjk2 * rik2;
        if (r6 > cutoff_triple_sixth) continue;

        nu_local = nu[type[i]][type[j]][type[k]];
        if (nu_local == 0.0) continue;

        interaction_ddd(nu_local, r6, rij2, rik2, rjk2,
                        rij, rik, rjk, fj, fk, EFLAG, evdwl);

        f[i].x -= fj[0] + fk[0];
        f[i].y -= fj[1] + fk[1];
        f[i].z -= fj[2] + fk[2];
        f[j].x += fj[0];
        f[j].y += fj[1];
        f[j].z += fj[2];
        f[k].x += fk[0];
        f[k].y += fk[1];
        f[k].z += fk[2];

        if (EVFLAG) ev_tally3_thr(this, i, j, k, evdwl, 0.0,
                                   fj, fk, rij, rik, thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairATMOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairATM::memory_usage();
  return bytes;
}

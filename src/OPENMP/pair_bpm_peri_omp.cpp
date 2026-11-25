/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

#include "pair_bpm_peri_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "lattice.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairBPMPeriOMP::PairBPMPeriOMP(LAMMPS *lmp) : PairBPMPeri(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairBPMPeriOMP::compute(int eflag, int vflag)
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
void PairBPMPeriOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int *_noalias const type = atom->type;
  const double *_noalias const special_lj = force->special_lj;
  const double *_noalias const vfrac = atom->dvector[index_vfrac];
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  const int nlocal = atom->nlocal;
  double evdwl = 0.0;

  // contact onset distance. Legacy uses min(0.9*r0_ref, 1.35*lc).
  // For non-bonded atoms (reference separation exceeds the horizon) this is always 1.35*lc,
  // => the dropped reference positions (x0) are not needed here
  const double d_ij = 1.35 * domain->lattice->xlattice;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    double fxtmp, fytmp, fztmp;
    fxtmp = fytmp = fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];

      // skip bonded pairs (excluded via the 1-2 special weight) and 1-3/1-4
      if (special_lj[sbmask(j)] == 0.0) continue;
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;
      const double r = sqrt(rsq);
      if (r >= d_ij) continue;

      // short-range repulsion (the 15 constant is from the EMU Theory Manual).
      // The mean nodal volume keeps the contact Newton-third-law balanced; for
      // uniform nodal volume this matches legacy's vfrac[j] form.
      const double dr = r - d_ij;
      const double vfrac_eff = 0.5 * (vfrac[i] + vfrac[j]);
      const double rk = 15.0 * kspring[itype][jtype] * vfrac_eff * (dr / cut[itype][jtype]);
      const double fpair = (r > 0.0) ? -(rk / r) : 0.0;

      fxtmp += delx * fpair;
      fytmp += dely * fpair;
      fztmp += delz * fpair;
      if (NEWTON_PAIR || (j < nlocal)) {
        f[j].x -= delx * fpair;
        f[j].y -= dely * fpair;
        f[j].z -= delz * fpair;
      }

      if (EFLAG) evdwl = 0.5 * rk * dr;
      // the peridynamic stress integrates over both nodal volumes, so the virial
      // carries one more vfrac factor than the (one-vfrac) contact force/energy
      // (matches legacy PERI's fpair*vfrac[i]); consistent with bond_style bpm/peri
      if (EVFLAG)
        ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, 0.0, fpair * vfrac_eff, delx, dely,
                     delz, thr);
    }

    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairBPMPeriOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBPMPeri::memory_usage();

  return bytes;
}

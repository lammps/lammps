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

#include "pair_coul_cut_dielectric_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;
using MathConst::MY_PIS;

static constexpr double EPSILON = 1.0e-6;

/* ---------------------------------------------------------------------- */

PairCoulCutDielectricOMP::PairCoulCutDielectricOMP(LAMMPS *_lmp) :
    PairCoulCutDielectric(_lmp), ThrOMP(_lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void PairCoulCutDielectricOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (atom->nmax > nmax) {
    memory->destroy(efield);
    nmax = atom->nmax;
    memory->create(efield, nmax, 3, "pair:efield");
  }

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
      if (eflag) eval<1,1>(ifrom, ito, thr);
      else       eval<1,0>(ifrom, ito, thr);
    } else {
      eval<0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG>
void PairCoulCutDielectricOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q_scaled;
  const double * _noalias const eps = atom->epsilon;
  const auto * _noalias const norm = (dbl4_t *) atom->mu[0];
  const double * _noalias const curvature = atom->curvature;
  const double * _noalias const area = atom->area;
  const int * _noalias const type = atom->type;
  const double * _noalias const special_coul = force->special_coul;
  const double qqrd2e = force->qqrd2e;

  const int * _noalias const ilist    = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  int ** const firstneigh             = list->firstneigh;

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double qtmp  = q[i];
    const double etmp  = eps[i];
    const double xtmp  = x[i].x;
    const double ytmp  = x[i].y;
    const double ztmp  = x[i].z;
    const int itype    = type[i];
    const double * _noalias const cutsqi  = cutsq[itype];
    const double * _noalias const scalei  = scale[itype];
    const int * _noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    double fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;
    double extmp = 0.0, eytmp = 0.0, eztmp = 0.0;

    // self term - thread-safe: outer atom i is owned by exactly one thread
    double curvature_threshold = sqrt(area[i]);
    if (curvature[i] < curvature_threshold) {
      double sf = curvature[i] / (4.0 * MY_PIS * curvature_threshold) * area[i] * q[i];
      efield[i][0] = sf * norm[i].x;
      efield[i][1] = sf * norm[i].y;
      efield[i][2] = sf * norm[i].z;
    } else {
      efield[i][0] = efield[i][1] = efield[i][2] = 0.0;
    }

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq  = delx*delx + dely*dely + delz*delz;
      const int jtype   = type[j];

      if (rsq < cutsqi[jtype] && rsq > EPSILON) {
        const double r2inv = 1.0/rsq;
        const double rinv  = sqrt(r2inv);

        double efield_i = qqrd2e * scalei[jtype] * q[j] * rinv;
        const double forcecoul = qtmp * efield_i;
        const double fpair_i = factor_coul * etmp * forcecoul * r2inv;

        fxtmp += delx * fpair_i;
        fytmp += dely * fpair_i;
        fztmp += delz * fpair_i;

        efield_i *= (factor_coul * etmp * r2inv);
        extmp += delx * efield_i;
        eytmp += dely * efield_i;
        eztmp += delz * efield_i;

        double ecoul = 0.0;
        if (EFLAG) {
          ecoul = factor_coul * qqrd2e * scalei[jtype] * qtmp * q[j] *
              0.5 * (etmp + eps[j]) * rinv;
        }

        if (EVFLAG) ev_tally_full_thr(this, i, 0.0, ecoul, fpair_i, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
    efield[i][0] += extmp;
    efield[i][1] += eytmp;
    efield[i][2] += eztmp;
  }
}

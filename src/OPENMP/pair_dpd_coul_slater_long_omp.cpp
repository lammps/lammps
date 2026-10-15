/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

#include "pair_dpd_coul_slater_long_omp.h"

#include "atom.h"
#include "comm.h"
#include "ewald_const.h"
#include "force.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "suffix.h"
#include "update.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace EwaldConst;

static constexpr double EPSILON = 1.0e-10;

/* ---------------------------------------------------------------------- */

PairDPDCoulSlaterLongOMP::PairDPDCoulSlaterLongOMP(LAMMPS *lmp) :
    PairDPDCoulSlaterLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  random_thr = nullptr;
  nthreads = 0;
}

/* ---------------------------------------------------------------------- */

PairDPDCoulSlaterLongOMP::~PairDPDCoulSlaterLongOMP()
{
  if (random_thr) {
    for (int i = 1; i < nthreads; ++i) delete random_thr[i];

    delete[] random_thr;
    random_thr = nullptr;
  }
}

/* ----------------------------------------------------------------------
   the base class settings() re-creates the serial random number generator,
   so the pool of per-thread generators must be re-created on the next call
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLongOMP::settings(int narg, char **arg)
{
  PairDPDCoulSlaterLong::settings(narg, arg);

  if (random_thr) {
    for (int i = 1; i < nthreads; ++i) delete random_thr[i];
    delete[] random_thr;
    random_thr = nullptr;
  }
  nthreads = 0;
}

/* ---------------------------------------------------------------------- */

void PairDPDCoulSlaterLongOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // precompute random force scaling factors

  for (int i = 0; i < 4; ++i) special_sqrt[i] = sqrt(force->special_lj[i]);

  const int nall = atom->nlocal + atom->nghost;
  const int inum = list->inum;

  // number of threads has changed. reallocate pool of pRNGs
  if (nthreads != comm->nthreads) {
    if (random_thr) {
      for (int i = 1; i < nthreads; ++i) delete random_thr[i];

      delete[] random_thr;
    }

    nthreads = comm->nthreads;
    random_thr = new RanMars *[nthreads];
    for (int i = 1; i < nthreads; ++i) random_thr[i] = nullptr;

    // to ensure full compatibility with the serial style
    // we use the serial random number generator instance for thread 0
    random_thr[0] = random;
  }

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    // generate a random number generator instance for
    // all threads != 0. make sure we use unique seeds.
    if ((tid > 0) && (random_thr[tid] == nullptr)) {
      // wrap around to stay within the valid seed range of RanMars
      const bigint tseed = (bigint) seed + comm->me + (bigint) comm->nprocs * tid;
      random_thr[tid] = new RanMars(Pair::lmp, (int) ((tseed - 1) % 900000000 + 1));
    }

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
void PairDPDCoulSlaterLongOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  const auto *_noalias const v = (dbl3_t *) atom->v[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const double *_noalias const q = atom->q;
  const int *_noalias const type = atom->type;
  const double *_noalias const special_lj = force->special_lj;
  const double *_noalias const special_coul = force->special_coul;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  const int nlocal = atom->nlocal;
  const double qqrd2e = force->qqrd2e;
  const double dtinvsqrt = 1.0 / sqrt(update->dt);
  RanMars &rng = *random_thr[thr->get_tid()];

  double evdwl = 0.0;
  double ecoul = 0.0;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    const double qtmp = q[i];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const double vxtmp = v[i].x;
    const double vytmp = v[i].y;
    const double vztmp = v[i].z;
    double fxtmp, fytmp, fztmp;
    fxtmp = fytmp = fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_dpd = special_lj[sbmask(j)];
      const double factor_sqrt = special_sqrt[sbmask(j)];
      const double factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];

      // forces if below maximum cutoff
      if (rsq < cutsq[itype][jtype]) {
        const double r = sqrt(rsq);
        if (EVFLAG) evdwl = ecoul = 0.0;
        double forcedpd, forcecoul;

        // apply DPD force if distance below DPD cutoff
        if ((rsq < cut_dpdsq[itype][jtype]) && (r > EPSILON)) {
          const double rinv = 1.0 / r;
          const double delvx = vxtmp - v[j].x;
          const double delvy = vytmp - v[j].y;
          const double delvz = vztmp - v[j].z;
          const double dot = delx * delvx + dely * delvy + delz * delvz;
          const double wd = 1.0 - r / cut_dpd[itype][jtype];
          const double randnum = rng.gaussian();

          // conservative force = a0 * wd
          // drag force = -gamma * wd^2 * (delx dot delv) / r
          // random force = sigma * wd * rnd * dtinvsqrt;
          // random force must be scaled by sqrt(factor_dpd)

          forcedpd = a0[itype][jtype] * wd;
          forcedpd -= gamma[itype][jtype] * wd * wd * dot * rinv;
          forcedpd *= factor_dpd;
          forcedpd += factor_sqrt * sigma[itype][jtype] * wd * randnum * dtinvsqrt;
          forcedpd *= rinv;

          if (EFLAG) {
            // eng shifted to 0.0 at cutoff
            evdwl = 0.5 * a0[itype][jtype] * cut_dpd[itype][jtype] * wd * wd;
            evdwl *= factor_dpd;
          }

        } else
          forcedpd = 0.0;

        // apply Slater electrostatic force if distance below Slater cutoff
        // and the two species are charged
        if (rsq < cut_slatersq[itype][jtype]) {
          const double r2inv = 1.0 / rsq;
          const double grij = g_ewald * r;
          const double expm2 = exp(-grij * grij);
          const double t = 1.0 / (1.0 + EWALD_P * grij);
          const double erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
          const double slater_term =
              exp(-2 * r / lamda) * (1 + (2 * r / lamda * (1 + r / lamda)));
          const double prefactor = qqrd2e * qtmp * q[j] / r;
          forcecoul = prefactor * (erfc + EWALD_F * grij * expm2 - slater_term);
          if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor * (1 - slater_term);
          forcecoul *= r2inv;

          if (EFLAG) {
            ecoul = prefactor * (erfc - (1 + r / lamda) * exp(-2 * r / lamda));
            if (factor_coul < 1.0)
              ecoul -= (1.0 - factor_coul) * prefactor *
                  (1.0 - (1 + r / lamda) * exp(-2 * r / lamda));
          }

        } else
          forcecoul = 0.0;

        const double fpair = forcedpd + forcecoul;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, ecoul, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairDPDCoulSlaterLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairDPDCoulSlaterLong::memory_usage();
  bytes += (double) comm->nthreads * sizeof(RanMars *);
  bytes += (double) comm->nthreads * sizeof(RanMars);

  return bytes;
}

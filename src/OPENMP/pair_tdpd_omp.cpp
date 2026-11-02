/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_tdpd_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "suffix.h"
#include "update.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

static constexpr double EPSILON = 1.0e-10;

/* ---------------------------------------------------------------------- */

PairTDPDOMP::PairTDPDOMP(LAMMPS *lmp) : PairTDPD(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  random_thr = nullptr;
  nthreads = 0;
}

/* ---------------------------------------------------------------------- */

PairTDPDOMP::~PairTDPDOMP()
{
  destroy_random_thr();
}

/* ----------------------------------------------------------------------
   thread 0 uses the generator of the base class
------------------------------------------------------------------------- */

void PairTDPDOMP::destroy_random_thr()
{
  if (random_thr) {
    for (int i = 1; i < nthreads; ++i) delete random_thr[i];
    delete[] random_thr;
    random_thr = nullptr;
  }
  nthreads = 0;
}

/* ----------------------------------------------------------------------
   the base class re-creates the serial random number generator here,
   so the pool of per-thread generators must be re-created on the next call
------------------------------------------------------------------------- */

void PairTDPDOMP::settings(int narg, char **arg)
{
  PairTDPD::settings(narg, arg);
  destroy_random_thr();
}

/* ---------------------------------------------------------------------- */

void PairTDPDOMP::read_restart_settings(FILE *fp)
{
  PairTDPD::read_restart_settings(fp);
  destroy_random_thr();
}

/* ---------------------------------------------------------------------- */

void PairTDPDOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int inum = list->inum;

  // number of threads has changed. reallocate pool of pRNGs
  if (nthreads != comm->nthreads) {
    destroy_random_thr();

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
void PairTDPDOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  const auto *_noalias const v = (dbl3_t *) atom->v[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  double *const *const cc_flux = thr->get_cc_flux();
  const double *const *const cc = atom->cc;
  const int *_noalias const type = atom->type;
  const double *_noalias const special_lj = force->special_lj;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;
  RanMars &rng = *random_thr[thr->get_tid()];

  const int nlocal = atom->nlocal;
  const double dtinvsqrt = 1.0 / sqrt(update->dt);
  double evdwl = 0.0;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const double vxtmp = v[i].x;
    const double vytmp = v[i].y;
    const double vztmp = v[i].z;
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    double fxtmp, fytmp, fztmp;
    fxtmp = fytmp = fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_dpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        const double r = sqrt(rsq);
        if (r < EPSILON) continue;    // r can be 0.0 in DPD systems
        const double rinv = 1.0 / r;
        const double delvx = vxtmp - v[j].x;
        const double delvy = vytmp - v[j].y;
        const double delvz = vztmp - v[j].z;
        const double dot = delx * delvx + dely * delvy + delz * delvz;
        double wc = 1.0 - r / cut[itype][jtype];
        wc = MAX(0, MIN(1.0, wc));
        const double wr = pow(wc, 0.5 * power[itype][jtype]);
        const double randnum = rng.gaussian();

        // conservative force = a0 * wc
        // drag force = -gamma * wr^2 * (delx dot delv) / r
        // random force = sigma * wr^(power/2) * rnd * dtinvsqrt;

        double fpair = a0[itype][jtype] * wc;
        fpair -= gamma[itype][jtype] * wr * wr * dot * rinv;
        fpair += sigma[itype][jtype] * wr * randnum * dtinvsqrt;
        fpair *= factor_dpd * rinv;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;

        // chemical concentration transport
        if (r < cutcc[itype][jtype]) {
          for (int k = 0; k < cc_species; k++) {
            double wcr = 1.0 - r / cutcc[itype][jtype];
            wcr = MAX(0, wcr);
            wcr = pow(wcr, 0.5 * powercc[itype][jtype][k]);
            double randnumcc = rng.gaussian();
            randnumcc = MAX(-5.0, MIN(randnumcc, 5.0));
            const double dQc = -kappa[itype][jtype][k] * wcr * wcr * (cc[i][k] - cc[j][k]);
            const double dQr =
                epsilon[itype][jtype][k] * wcr * (cc[i][k] + cc[j][k]) * randnumcc * dtinvsqrt;
            cc_flux[i][k] += (dQc + dQr);
            if (NEWTON_PAIR || j < nlocal) cc_flux[j][k] -= (dQc + dQr);
          }
        }

        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        if (EFLAG) {
          evdwl = 0.5 * a0[itype][jtype] * cut[itype][jtype] * wc * wc;
          evdwl *= factor_dpd;
        }

        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, 0.0, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairTDPDOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairTDPD::memory_usage();

  return bytes;
}

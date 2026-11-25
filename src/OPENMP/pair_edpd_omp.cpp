/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_edpd_omp.h"

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

PairEDPDOMP::PairEDPDOMP(LAMMPS *lmp) : PairEDPD(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  random_thr = randomT_thr = nullptr;
  nthreads = 0;
}

/* ---------------------------------------------------------------------- */

PairEDPDOMP::~PairEDPDOMP()
{
  destroy_random_thr();
}

/* ----------------------------------------------------------------------
   thread 0 uses the generators of the base class
------------------------------------------------------------------------- */

void PairEDPDOMP::destroy_random_thr()
{
  if (random_thr) {
    for (int i = 1; i < nthreads; ++i) {
      delete random_thr[i];
      delete randomT_thr[i];
    }
    delete[] random_thr;
    delete[] randomT_thr;
    random_thr = randomT_thr = nullptr;
  }
  nthreads = 0;
}

/* ----------------------------------------------------------------------
   the base class re-creates the serial random number generators here,
   so the pools of per-thread generators must be re-created on the next call
------------------------------------------------------------------------- */

void PairEDPDOMP::settings(int narg, char **arg)
{
  PairEDPD::settings(narg, arg);
  destroy_random_thr();
}

/* ---------------------------------------------------------------------- */

void PairEDPDOMP::read_restart_settings(FILE *fp)
{
  PairEDPD::read_restart_settings(fp);
  destroy_random_thr();
}

/* ---------------------------------------------------------------------- */

void PairEDPDOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int inum = list->inum;

  // number of threads has changed. reallocate pools of pRNGs
  if (nthreads != comm->nthreads) {
    destroy_random_thr();

    nthreads = comm->nthreads;
    random_thr = new RanMars *[nthreads];
    randomT_thr = new RanMars *[nthreads];
    for (int i = 1; i < nthreads; ++i) random_thr[i] = randomT_thr[i] = nullptr;

    // to ensure full compatibility with the serial style
    // we use the serial random number generator instances for thread 0
    random_thr[0] = random;
    randomT_thr[0] = randomT;
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

    // generate random number generator instances for
    // all threads != 0. make sure we use unique seeds.
    if ((tid > 0) && (random_thr[tid] == nullptr)) {
      // wrap around to stay within the valid seed range of RanMars
      const bigint offset = comm->me + (bigint) comm->nprocs * tid;
      random_thr[tid] =
          new RanMars(Pair::lmp, (int) (((bigint) seed + offset - 1) % 900000000 + 1));
      randomT_thr[tid] =
          new RanMars(Pair::lmp, (int) ((2 * (bigint) seed + offset - 1) % 900000000 + 1));
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
void PairEDPDOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  const auto *_noalias const v = (dbl3_t *) atom->v[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  double *_noalias const Q = thr->get_edpd_flux();
  const double *_noalias const T = atom->edpd_temp;
  const double *_noalias const cv = atom->edpd_cv;
  const double *_noalias const mass = atom->mass;
  const int *_noalias const type = atom->type;
  const double *_noalias const special_lj = force->special_lj;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;
  RanMars &rng = *random_thr[thr->get_tid()];
  RanMars &rngT = *randomT_thr[thr->get_tid()];

  const int nlocal = atom->nlocal;
  const double dtinvsqrt = 1.0 / sqrt(update->dt);
  const double kboltz = 1.0;
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
        if (r < EPSILON) continue;
        const double rinv = 1.0 / r;
        const double delvx = vxtmp - v[j].x;
        const double delvy = vytmp - v[j].y;
        const double delvz = vztmp - v[j].z;
        const double dot = delx * delvx + dely * delvy + delz * delvz;
        const double vijeij = dot * rinv;
        const double randnum = rng.gaussian();

        const double T_ij = 0.5 * (T[i] + T[j]);
        double T_pow[4];
        T_pow[0] = T_ij - 1.0;
        T_pow[1] = T_pow[0] * T_pow[0];
        T_pow[2] = T_pow[0] * T_pow[1];
        T_pow[3] = T_pow[0] * T_pow[2];

        double power_d = power[itype][jtype];
        if (power_flag) {
          double factor = 1.0;
          for (int k = 0; k < 4; k++) factor += sc[itype][jtype][k] * T_pow[k];
          power_d *= factor;
        }

        power_d = MAX(0.01, power_d);
        double wc = 1.0 - r / cut[itype][jtype];
        wc = MAX(0.0, MIN(1.0, wc));
        const double wr = pow(wc, 0.5 * power_d);

        const double GammaIJ = gamma[itype][jtype];
        const double SigmaIJ = sqrt(4.0 * GammaIJ * kboltz * T[i] * T[j] / (T[i] + T[j]));

        double fpair = a0[itype][jtype] * T_ij * wc;
        fpair -= GammaIJ * wr * wr * dot * rinv;
        fpair += SigmaIJ * wr * randnum * dtinvsqrt;
        fpair *= factor_dpd * rinv;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;

        // heat transfer
        double dQc = 0.0, dQd = 0.0, dQr = 0.0;
        if (r < cutT[itype][jtype]) {
          double wrT = 1.0 - r / cutT[itype][jtype];
          wrT = MAX(0.0, MIN(1.0, wrT));
          wrT = pow(wrT, 0.5 * powerT[itype][jtype]);
          double randnumT = rngT.gaussian();
          randnumT = MAX(-5.0, MIN(randnumT, 5.0));

          double kappaT = kappa[itype][jtype];
          if (kappa_flag) {
            double factor = 1.0;
            for (int k = 0; k < 4; k++) factor += kc[itype][jtype][k] * T_pow[k];
            kappaT *= factor;
          }

          const double kij = cv[i] * cv[j] * kappaT * T_ij * T_ij;
          const double alphaij = sqrt(2.0 * kboltz * kij);

          dQc = kij * wrT * wrT * (T[j] - T[i]) / (T[i] * T[j]);
          dQd = wr * wr * (GammaIJ * vijeij * vijeij - SigmaIJ * SigmaIJ / mass[itype]) -
              SigmaIJ * wr * vijeij * randnum;
          dQd /= (cv[i] + cv[j]);
          dQr = alphaij * wrT * dtinvsqrt * randnumT;
          Q[i] += (dQc + dQd + dQr);
        }

        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
          Q[j] -= (dQc - dQd + dQr);
        }

        if (EFLAG) {
          evdwl = 0.5 * a0[itype][jtype] * T_ij * cut[itype][jtype] * wc * wc;
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

double PairEDPDOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairEDPD::memory_usage();

  return bytes;
}

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_sdpd_taitwater_isothermal_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "suffix.h"
#include "update.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

static const double sqrt_2_inv = std::sqrt(0.5);

/* ---------------------------------------------------------------------- */

PairSDPDTaitwaterIsothermalOMP::PairSDPDTaitwaterIsothermalOMP(LAMMPS *lmp) :
    PairSDPDTaitwaterIsothermal(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  random_thr = nullptr;
  nthreads = 0;
}

/* ---------------------------------------------------------------------- */

PairSDPDTaitwaterIsothermalOMP::~PairSDPDTaitwaterIsothermalOMP()
{
  destroy_random_thr();
}

/* ----------------------------------------------------------------------
   thread 0 uses the generator of the base class, unless that is the
   (optional) zest generator, which the threads cannot share
------------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermalOMP::destroy_random_thr()
{
  if (random_thr) {
#ifdef USE_ZEST
    const int first = 0;
#else
    const int first = 1;
#endif
    for (int i = first; i < nthreads; ++i) delete random_thr[i];
    delete[] random_thr;
    random_thr = nullptr;
  }
  nthreads = 0;
}

/* ----------------------------------------------------------------------
   the base class settings() re-creates the serial random number generator,
   so the pool of per-thread generators must be re-created on the next call
------------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermalOMP::settings(int narg, char **arg)
{
  PairSDPDTaitwaterIsothermal::settings(narg, arg);
  destroy_random_thr();
}

/* ---------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermalOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int inum = list->inum;

  // number of threads has changed. reallocate pool of pRNGs
  if (nthreads != comm->nthreads) {
    destroy_random_thr();

    nthreads = comm->nthreads;
    random_thr = new RanMars *[nthreads];
    for (int i = 0; i < nthreads; ++i) random_thr[i] = nullptr;

    // to ensure full compatibility with the serial style
    // we use the serial random number generator instance for thread 0
#ifndef USE_ZEST
    random_thr[0] = random;
#endif
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
    if (random_thr[tid] == nullptr) {
      // wrap around to stay within the valid seed range of RanMars
      const bigint tseed = (bigint) seed + (bigint) comm->nprocs * tid;
      random_thr[tid] = new RanMars(Pair::lmp, (int) ((tseed - 1) % 900000000 + 1));
    }

    // this style tallies no energy

    if (evflag) {
      if (force->newton_pair)
        eval<1, 1>(ifrom, ito, thr);
      else
        eval<1, 0>(ifrom, ito, thr);
    } else {
      if (force->newton_pair)
        eval<0, 1>(ifrom, ito, thr);
      else
        eval<0, 0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int NEWTON_PAIR>
void PairSDPDTaitwaterIsothermalOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  const auto *_noalias const v = (dbl3_t *) atom->vest[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  double *_noalias const drho = thr->get_drho();
  const double *_noalias const rho = atom->rho;
  const double *_noalias const mass = atom->mass;
  const int *_noalias const type = atom->type;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;
  RanMars &rng = *random_thr[thr->get_tid()];

  const int nlocal = atom->nlocal;
  const int dimension = domain->dimension;
  const double dtinv = 1.0 / update->dt;
  const double kBoltzmann = force->boltz;
  double wiener[3][3], f_random[3];

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

    const double imass = mass[itype];

    // compute pressure of atom i with Tait EOS
    double tmp = rho[i] / rho0[itype];
    double fi = tmp * tmp * tmp;
    fi = B[itype] * (fi * fi * tmp - 1.0) / (rho[i] * rho[i]);

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];
      const double jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
        const double h = cut[itype][jtype];
        const double ih = 1.0 / h;
        const double ihsq = ih * ih;

        const double r = sqrt(rsq);
        double wfd = h - r;
        if (dimension == 3) {
          // Lucy Kernel, 3d
          // wfd lacks a factor of r, see the base style
          wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        }

        // compute pressure of atom j with Tait EOS
        tmp = rho[j] / rho0[jtype];
        double fj = tmp * tmp * tmp;
        fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[j] * rho[j]);

        const double velx = vxtmp - v[j].x;
        const double vely = vytmp - v[j].y;
        const double velz = vztmp - v[j].z;

        // dot product of velocity delta and distance vector
        const double delVdotDelR = delx * velx + dely * vely + delz * velz;

        // Espanol Viscosity (Espanol, 2003)

        const double fvisc = (5. / 3.) * viscosity * imass * jmass * wfd / (rho[i] * rho[j]);

        // total pair force
        const double fpair = -imass * jmass * (fi + fj) * wfd;

        // random force calculation
        // independent increments of a Wiener process matrix
        wiener[0][0] = rng.gaussian();
        wiener[1][1] = rng.gaussian();
        wiener[2][2] = rng.gaussian();

        wiener[0][1] = wiener[1][0] = sqrt_2_inv * rng.gaussian();
        wiener[0][2] = wiener[2][0] = sqrt_2_inv * rng.gaussian();
        wiener[1][2] = wiener[2][1] = sqrt_2_inv * rng.gaussian();

        const double prefactor = sqrt(-4. * kBoltzmann * temperature * fvisc * dtinv) / r;

        f_random[0] = prefactor * (wiener[0][0] * delx + wiener[0][1] * dely + wiener[0][2] * delz);
        f_random[1] = prefactor * (wiener[1][0] * delx + wiener[1][1] * dely + wiener[1][2] * delz);
        f_random[2] = prefactor * (wiener[2][0] * delx + wiener[2][1] * dely + wiener[2][2] * delz);

        f[i].x += delx * fpair + (velx + delx * delVdotDelR / rsq) * fvisc + f_random[0];
        f[i].y += dely * fpair + (vely + dely * delVdotDelR / rsq) * fvisc + f_random[1];
        f[i].z += delz * fpair + (velz + delz * delVdotDelR / rsq) * fvisc + f_random[2];

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair + (velx + delx * delVdotDelR / rsq) * fvisc + f_random[0];
          f[j].y -= dely * fpair + (vely + dely * delVdotDelR / rsq) * fvisc + f_random[1];
          f[j].z -= delz * fpair + (velz + delz * delVdotDelR / rsq) * fvisc + f_random[2];
          drho[j] += imass * delVdotDelR * wfd;
        }

        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, 0.0, 0.0, fpair, delx, dely, delz, thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairSDPDTaitwaterIsothermalOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairSDPDTaitwaterIsothermal::memory_usage();

  return bytes;
}

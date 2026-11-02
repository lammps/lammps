/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

#include "pair_sph_taitwater_morris_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHTaitwaterMorrisOMP::PairSPHTaitwaterMorrisOMP(LAMMPS *lmp) :
    PairSPHTaitwaterMorris(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairSPHTaitwaterMorrisOMP::compute(int eflag, int vflag)
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
void PairSPHTaitwaterMorrisOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  const auto *_noalias const v = (dbl3_t *) atom->vest[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  double *_noalias const desph = thr->get_de();
  double *_noalias const drho = thr->get_drho();
  const double *_noalias const rho = atom->rho;
  const double *_noalias const mass = atom->mass;
  const int *_noalias const type = atom->type;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  const int nlocal = atom->nlocal;
  const int dimension = domain->dimension;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const double vxtmp = v[i].x;
    const double vytmp = v[i].y;
    const double vztmp = v[i].z;
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

        double wfd = h - sqrt(rsq);
        if (dimension == 3) {
          // Lucy Kernel, 3d
          // Note that wfd, the derivative of the weight function with respect to r,
          // is lacking a factor of r.
          // The missing factor of r is recovered by
          // (1) using delV . delX instead of delV . (delX/r) and
          // (2) using f[i][0] += delx * fpair instead of f[i][0] += (delx/r) * fpair
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

        // Morris Viscosity (Morris, 1996)

        double fvisc = 2 * viscosity[itype][jtype] / (rho[i] * rho[j]);

        fvisc *= imass * jmass * wfd;

        // total pair force & thermal energy increment
        const double fpair = -imass * jmass * (fi + fj) * wfd;
        const double deltaE =
            -0.5 * (fpair * delVdotDelR + fvisc * (velx * velx + vely * vely + velz * velz));

        f[i].x += delx * fpair + velx * fvisc;
        f[i].y += dely * fpair + vely * fvisc;
        f[i].z += delz * fpair + velz * fvisc;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        // change in thermal energy
        desph[i] += deltaE;

        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair + velx * fvisc;
          f[j].y -= dely * fpair + vely * fvisc;
          f[j].z -= delz * fpair + velz * fvisc;
          desph[j] += deltaE;
          drho[j] += imass * delVdotDelR * wfd;
        }

        // viscous forces do not contribute to virial
        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, 0.0, 0.0, fpair, delx, dely, delz, thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairSPHTaitwaterMorrisOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairSPHTaitwaterMorris::memory_usage();

  return bytes;
}

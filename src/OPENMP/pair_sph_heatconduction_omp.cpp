/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Claude Fable 5.1 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_sph_heatconduction_omp.h"

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

PairSPHHeatConductionOMP::PairSPHHeatConductionOMP(LAMMPS *lmp) :
    PairSPHHeatConduction(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairSPHHeatConductionOMP::compute(int eflag, int vflag)
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

    if (force->newton_pair)
      eval<1>(ifrom, ito, thr);
    else
      eval<0>(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int NEWTON_PAIR>
void PairSPHHeatConductionOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  double *_noalias const desph = thr->get_de();
  const double *_noalias const esph = atom->esph;
  const double *_noalias const rho = atom->rho;
  const double *_noalias const mass = atom->mass;
  const int *_noalias const type = atom->type;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  const int nlocal = atom->nlocal;
  const int dimension = domain->dimension;

  // loop over neighbors of my atoms and do heat diffusion

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const double imass = mass[itype];

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

        // kernel function
        double wfd = h - sqrt(rsq);
        if (dimension == 3) {
          // Lucy Kernel, 3d
          // Note that wfd, the derivative of the weight function with respect to r,
          // is lacking a factor of r.
          // The missing factor of r is recovered by
          // deltaE, which is missing a factor of 1/r
          wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        }

        const double D = alpha[itype][jtype];    // diffusion coefficient

        double deltaE = 2.0 * imass * jmass / (imass + jmass);
        deltaE *= (rho[i] + rho[j]) / (rho[i] * rho[j]);
        deltaE *= D * (esph[i] - esph[j]) * wfd;

        desph[i] += deltaE;
        if (NEWTON_PAIR || j < nlocal) desph[j] -= deltaE;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairSPHHeatConductionOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairSPHHeatConduction::memory_usage();

  return bytes;
}

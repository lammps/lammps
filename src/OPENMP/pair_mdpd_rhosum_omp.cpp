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

#include "pair_mdpd_rhosum_omp.h"

#include "atom.h"
#include "comm.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMDPDRhoSumOMP::PairMDPDRhoSumOMP(LAMMPS *lmp) : PairMDPDRhoSum(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairMDPDRhoSumOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // recompute density using the full neighbor list.
  // each thread only writes the density of its own atoms, so no reduction is needed.

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    eval(ifrom, ito);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region

  // communicate densities
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void PairMDPDRhoSumOMP::eval(int iifrom, int iito)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  double *_noalias const rho = atom->rho;
  const double *_noalias const mass = atom->mass;
  const int *_noalias const type = atom->type;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;

    // initialize density with zero self-contribution

    rho[i] = 0.0;

    // add density at each atom via kernel function overlap

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      const int jtype = type[j];
      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[itype][jtype]) {
        const double h = cut[itype][jtype];
        const double ih = 1.0 / h;
        const double ihsq = ih * ih;

        // Lucy kernel, 3d
        const double r = sqrt(rsq);
        double wf = (h - r) * ihsq;
        wf = 2.0889086280811262819e0 * (h + 3. * r) * wf * wf * wf * ih;
        rho[i] += mass[jtype] * wf;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairMDPDRhoSumOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairMDPDRhoSum::memory_usage();

  return bytes;
}

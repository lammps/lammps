/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

#include "pair_multi_lucy_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

PairMultiLucyOMP::PairMultiLucyOMP(LAMMPS *lmp) :
    PairMultiLucy(lmp), ThrOMP(lmp, THR_PAIR), rho_thr(nullptr), nmax_thr(0)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairMultiLucyOMP::~PairMultiLucyOMP()
{
  memory->destroy(rho_thr);
}

/* ---------------------------------------------------------------------- */

void PairMultiLucyOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow the per-thread copies of the local density array if necessary

  if (atom->nmax > nmax_thr) {
    nmax_thr = atom->nmax;
    memory->grow(rho_thr, nmax_thr * nthreads, "pair_multi_lucy:rho_thr");
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

    // hook up and zero the per-thread copy of the local density array
    thr->init_eam(nall, rho_thr);

    if (evflag) {
      if (eflag) {
        eval<1, 1>(ifrom, ito, thr);
      } else {
        eval<1, 0>(ifrom, ito, thr);
      }
    } else {
      eval<0, 0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region

  error_thr();
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG>
void PairMultiLucyOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  double *const rho_t = thr->get_rho();
  double *const rho = atom->rho;
  const int *_noalias const type = atom->type;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  const int tid = thr->get_tid();
  const int nthreads = comm->nthreads;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const int newton_pair = force->newton_pair;
  const int tlm1 = tablength - 1;

  // local density at each atom from the Lucy kernel, accumulated per thread
  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      const int jtype = type[j];

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[itype][jtype]) {
        const double rcut = sqrt(cutsq[itype][jtype]);
        const double factor = (84.0 / (5.0 * MY_PI * rcut * rcut * rcut)) *
            (1.0 + 3.0 * sqrt(rsq) / (2.0 * rcut)) * (1.0 - sqrt(rsq) / rcut) *
            (1.0 - sqrt(rsq) / rcut) * (1.0 - sqrt(rsq) / rcut) * (1.0 - sqrt(rsq) / rcut);
        rho_t[i] += factor;
        if (newton_pair || j < nlocal) rho_t[j] += factor;
      }
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // reduce the per-thread copies of the local density
  thr->timer(Timer::PAIR);
  data_reduce_thr(rho_thr, nall, nthreads, 1, tid);

  // wait until reduction is complete
  sync_threads();

  // copy the reduced density to the per-atom array and communicate it
  // MPI communication only on master thread
#if defined(_OPENMP)
#if _OPENMP >= 202011
#pragma omp masked
#else
#pragma omp master
#endif
#endif
  {
    memcpy(rho, rho_thr, nall * sizeof(double));
    if (newton_pair) comm->reverse_comm(this);
    comm->forward_comm(this);
  }

  // wait until master thread is done with communication
  sync_threads();

  // pairwise forces from the density dependent tables and
  // the density dependent self energy of each atom

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    const double rhoisq = rho[i] * rho[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        const Table *tb = &tables[tabindex[itype][jtype]];
        const double rhojsq = rho[j] * rho[j];
        if (check_error_thr(((rhoisq < tb->innersq) || (rhojsq < tb->innersq)), tid, FLERR,
                            "Density < table inner cutoff"))
          return;

        double fpair, A_i, A_j;
        if (tabstyle == LOOKUP) {
          const int itable = static_cast<int>((rhoisq - tb->innersq) * tb->invdelta);
          const int jtable = static_cast<int>((rhojsq - tb->innersq) * tb->invdelta);
          if (check_error_thr(((itable >= tlm1) || (jtable >= tlm1)), tid, FLERR,
                              "Density > table outer cutoff"))
            return;

          A_i = tb->f[itable];
          A_j = tb->f[jtable];
          fpair = 0.5 * (A_i + A_j) * (1.0 + 3.0 * sqrt(rsq) / sqrt(cutsq[itype][jtype])) *
              (1.0 - sqrt(rsq) / sqrt(cutsq[itype][jtype])) *
              (1.0 - sqrt(rsq) / sqrt(cutsq[itype][jtype])) *
              (1.0 - sqrt(rsq) / sqrt(cutsq[itype][jtype]));
          fpair = fpair / sqrt(rsq);

        } else if (tabstyle == LINEAR) {
          const int itable = static_cast<int>((rhoisq - tb->innersq) * tb->invdelta);
          const int jtable = static_cast<int>((rhojsq - tb->innersq) * tb->invdelta);
          if (check_error_thr(((itable >= tlm1) || (jtable >= tlm1)), tid, FLERR,
                              "Density > table outer cutoff"))
            return;

          const double fraction_i = ((rhoisq - tb->rsq[itable]) * tb->invdelta);
          const double fraction_j = ((rhojsq - tb->rsq[jtable]) * tb->invdelta);

          A_i = tb->f[itable] + fraction_i * tb->df[itable];
          A_j = tb->f[jtable] + fraction_j * tb->df[jtable];

          fpair = 0.5 * (A_i + A_j) * (1.0 + 3.0 * sqrt(rsq) / sqrt(cutsq[itype][jtype])) *
              (1.0 - sqrt(rsq) / sqrt(cutsq[itype][jtype])) *
              (1.0 - sqrt(rsq) / sqrt(cutsq[itype][jtype])) *
              (1.0 - sqrt(rsq) / sqrt(cutsq[itype][jtype]));

          fpair = fpair / sqrt(rsq);

        } else {
          if (check_error_thr(true, tid, FLERR,
                              "Only LOOKUP and LINEAR table styles have been implemented "
                              "for pair multi/lucy"))
            return;
          fpair = 0.0;
        }

        f[i].x += delx * fpair;
        f[i].y += dely * fpair;
        f[i].z += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }
        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz, thr);
      }
    }

    const Table *tb = &tables[tabindex[itype][itype]];
    if (check_error_thr((rhoisq < tb->innersq), tid, FLERR, "Density < table inner cutoff"))
      return;
    const int itable = static_cast<int>((rhoisq - tb->innersq) * tb->invdelta);
    double evdwl;
    if (tabstyle == LOOKUP)
      evdwl = tb->e[itable];
    else if (tabstyle == LINEAR) {
      if (check_error_thr((itable >= tlm1), tid, FLERR, "Density > table outer cutoff")) return;
      double fraction_i;
      if (itable == 0)
        fraction_i = 0.0;
      else
        fraction_i = ((rhoisq - tb->rsq[itable]) * tb->invdelta);
      evdwl = tb->e[itable] + fraction_i * tb->de[itable];
    } else {
      if (check_error_thr(true, tid, FLERR,
                          "Only LOOKUP and LINEAR table styles have been implemented "
                          "for pair multi/lucy"))
        return;
      evdwl = 0.0;
    }

    evdwl *= (MY_PI * cutsq[itype][itype] * cutsq[itype][itype]) / 84.0;

    if (EVFLAG) e_tally_thr(this, i, i, nlocal, newton_pair, evdwl, 0.0, thr);
  }
}

/* ---------------------------------------------------------------------- */

double PairMultiLucyOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairMultiLucy::memory_usage();
  bytes += (double) nmax_thr * comm->nthreads * sizeof(double);

  return bytes;
}

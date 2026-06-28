// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_qeq_dynamic_omp.h"

#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "neigh_list.h"

#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "omp_compat.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixQEqDynamicOMP::FixQEqDynamicOMP(LAMMPS *lmp, int narg, char **arg) :
  FixQEqDynamic(lmp, narg, arg), b_temp(nullptr), nmax_btmp(0)
{
}

/* ---------------------------------------------------------------------- */

FixQEqDynamicOMP::~FixQEqDynamicOMP()
{
  memory->destroy(b_temp);
}

/* ----------------------------------------------------------------------
   threaded matrix-free charge-force evaluation.  The forward/reverse
   communication (which zeroes ghost qf and sums ghost contributions back
   to their owners) stays serial; only the O(N*neigh) accumulation loop is
   threaded.  The qf[j] cross-write is deferred to per-thread b_temp and
   reduced, reproducing the serial FixQEqDynamic::compute_eneg().
------------------------------------------------------------------------- */

double FixQEqDynamicOMP::compute_eneg()
{
  int *ilist, *numneigh, **firstneigh;
  int *type = atom->type;
  int *mask = atom->mask;
  double *q = atom->q;
  double **x = atom->x;

  const int inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  const int nthreads = comm->nthreads;
  if (atom->nmax > nmax_btmp) {
    memory->destroy(b_temp);
    nmax_btmp = atom->nmax;
    memory->create(b_temp, nthreads, nmax_btmp, "qeq/dynamic/omp:b_temp");
  }

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) qf[i] = 0.0;
  }

  // communicate charge force to all nodes, first forward then reverse
  pack_flag = 2;
  comm->forward_comm(this);

  const int nall = atom->nlocal + atom->nghost;

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int tid = 0;
#if defined(_OPENMP)
    tid = omp_get_thread_num();
#endif
    double *qf_t = b_temp[tid];
    for (int i = 0; i < nall; ++i) qf_t[i] = 0.0;

#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50)
#endif
    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const int itype = type[i];
      if (mask[i] & groupbit) {
        qf[i] += chi[itype] + eta[itype] * q[i];

        const int *jlist = firstneigh[i];
        const int jnum = numneigh[i];

        for (int jj = 0; jj < jnum; jj++) {
          const int j = jlist[jj] & NEIGHMASK;
          const double dx = x[i][0] - x[j][0];
          const double dy = x[i][1] - x[j][1];
          const double dz = x[i][2] - x[j][2];
          const double rsq = dx*dx + dy*dy + dz*dz;
          if (rsq > cutoff_sq) continue;

          const double rinv = 1.0 / sqrt(rsq);
          qf[i] += q[j] * rinv;
          qf_t[j] += q[i] * rinv;
        }
      }
    }

#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50)
#endif
    for (int i = 0; i < nall; ++i)
      for (int t = 0; t < nthreads; ++t) qf[i] += b_temp[t][i];
  }

  pack_flag = 2;
  comm->reverse_comm(this);

  // sum charge force on each node and return it

  double eneg = 0.0, enegtot = 0.0;
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) eneg += qf[i];
  }
  MPI_Allreduce(&eneg, &enegtot, 1, MPI_DOUBLE, MPI_SUM, world);
  return enegtot;
}

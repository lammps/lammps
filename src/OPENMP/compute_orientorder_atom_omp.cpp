/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_orientorder_atom_omp.h"

#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeOrientOrderAtomOMP::ComputeOrientOrderAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeOrientOrderAtom(lmp, narg, arg), nthreads_alloc(0), maxneigh_thr(0), distsq_thr(nullptr),
    nearest_thr(nullptr), rlist_thr(nullptr), qnm_r_thr(nullptr), qnm_i_thr(nullptr)
{
}

/* ---------------------------------------------------------------------- */

ComputeOrientOrderAtomOMP::~ComputeOrientOrderAtomOMP()
{
  free_thr_arrays();
}

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtomOMP::free_thr_arrays()
{
  memory->destroy(distsq_thr);
  memory->destroy(nearest_thr);
  memory->destroy(rlist_thr);
  memory->destroy(qnm_r_thr);
  memory->destroy(qnm_i_thr);
  distsq_thr = nullptr;
  nearest_thr = nullptr;
  rlist_thr = nullptr;
  qnm_r_thr = nullptr;
  qnm_i_thr = nullptr;
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeOrientOrderAtom::compute_peratom().  The
   per-atom scratch arrays (distsq/nearest/rlist) and the boop accumulators
   (qnm_r/qnm_i) are replicated per thread and sized once, before the
   threaded region, to the largest neighbor count over the looped atoms.
   The base select3()/calc_boop() helpers are reentrant given thread-private
   scratch, so each atom is processed independently and the result is
   bit-identical to the serial compute.
------------------------------------------------------------------------- */

void ComputeOrientOrderAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow order parameter array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(qnarray);
    nmax = atom->nmax;
    memory->create(qnarray, nmax, ncol, "orientorder/atom:qnarray");
    array_atom = qnarray;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  double **x = atom->x;
  const int *const mask = atom->mask;
  memset(&qnarray[0][0], 0, sizeof(double) * nmax * ncol);

  // (re)allocate per-thread boop accumulators when the thread count changes

  const int nthreads = comm->nthreads;
  if (nthreads != nthreads_alloc) {
    free_thr_arrays();
    nthreads_alloc = nthreads;
    maxneigh_thr = 0;
    memory->create(qnm_r_thr, nthreads_alloc, nqlist, qmax + 1, "orientorder/atom:qnm_r_thr");
    memory->create(qnm_i_thr, nthreads_alloc, nqlist, qmax + 1, "orientorder/atom:qnm_i_thr");
  }

  // size the per-thread neighbor scratch once, to the largest neighbor count
  // (clamp to >= 1 so the buffers are never null when the region runs)

  int maxn = 1;
  for (int ii = 0; ii < inum; ii++) {
    const int jnum = numneigh[ilist[ii]];
    if (jnum > maxn) maxn = jnum;
  }
  if (maxn > maxneigh_thr) {
    memory->destroy(distsq_thr);
    memory->destroy(nearest_thr);
    memory->destroy(rlist_thr);
    maxneigh_thr = maxn;
    memory->create(distsq_thr, nthreads_alloc, maxneigh_thr, "orientorder/atom:distsq_thr");
    memory->create(nearest_thr, nthreads_alloc, maxneigh_thr, "orientorder/atom:nearest_thr");
    memory->create(rlist_thr, nthreads_alloc, maxneigh_thr, 3, "orientorder/atom:rlist_thr");
  }

  // compute order parameter for each atom in group
  // use full neighbor list to count atoms less than cutoff

#if defined(_OPENMP)
#pragma omp parallel
#endif
  {
    int tid = 0;
#if defined(_OPENMP)
    tid = omp_get_thread_num();
#endif
    double *const t_distsq = distsq_thr[tid];
    int *const t_nearest = nearest_thr[tid];
    double **const t_rlist = rlist_thr[tid];
    double **const t_qnm_r = qnm_r_thr[tid];
    double **const t_qnm_i = qnm_i_thr[tid];

#if defined(_OPENMP)
#pragma omp for schedule(dynamic)
#endif
    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      double *const qn = qnarray[i];
      if (mask[i] & groupbit) {
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        const int *const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        // loop over list of all neighbors within force cutoff
        // t_distsq[] = distance sq to each
        // t_rlist[] = distance vector to each
        // t_nearest[] = atom indices of neighbors

        int ncount = 0;
        for (int jj = 0; jj < jnum; jj++) {
          const int j = jlist[jj] & NEIGHMASK;

          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            t_distsq[ncount] = rsq;
            t_rlist[ncount][0] = delx;
            t_rlist[ncount][1] = dely;
            t_rlist[ncount][2] = delz;
            t_nearest[ncount++] = j;
          }
        }

        // if not nnn neighbors, order parameter = 0;

        if ((ncount == 0) || (ncount < nnn)) {
          for (int jj = 0; jj < ncol; jj++) qn[jj] = 0.0;
          continue;
        }

        // if nnn > 0, use only nearest nnn neighbors

        if (nnn > 0) {
          select3(nnn, ncount, t_distsq, t_nearest, t_rlist);
          ncount = nnn;
        }

        calc_boop(t_rlist, ncount, qn, qlist, nqlist, t_qnm_r, t_qnm_i);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array plus per-thread scratch
------------------------------------------------------------------------- */

double ComputeOrientOrderAtomOMP::memory_usage()
{
  double bytes = ComputeOrientOrderAtom::memory_usage();
  if (nthreads_alloc > 0) {
    bytes += (double) nthreads_alloc * maxneigh_thr * (4 * sizeof(double) + sizeof(int));
    bytes += (double) nthreads_alloc * nqlist * (qmax + 1) * 2 * sizeof(double);
  }
  return bytes;
}

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

#include "compute_hexorder_atom_omp.h"

#include "atom.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeHexOrderAtomOMP::ComputeHexOrderAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeHexOrderAtom(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeHexOrderAtom::compute_peratom().  The per-atom
   neighbor scratch (distsq/nearest) is allocated thread-locally inside the
   parallel region, sized to the largest neighbor count over the looped
   atoms.  select2()/calc_qn_complex() are reentrant given thread-private
   scratch, and each atom writes only its own output row qnarray[i], so the
   result is bit-identical to the serial compute.
------------------------------------------------------------------------- */

void ComputeHexOrderAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow order parameter array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(qnarray);
    nmax = atom->nmax;
    memory->create(qnarray, nmax, ncol, "hexorder/atom:qnarray");
    array_atom = qnarray;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // compute order parameter for each atom in group
  // use full neighbor list to count atoms less than cutoff

  const double *const *const x = atom->x;
  const int *const mask = atom->mask;

  // largest neighbor count over the looped atoms, for sizing thread scratch

  int maxn = 1;
  for (int ii = 0; ii < inum; ii++) {
    const int jnum = numneigh[ilist[ii]];
    if (jnum > maxn) maxn = jnum;
  }

#if defined(_OPENMP)
#pragma omp parallel
#endif
  {
    auto *const t_distsq = new double[maxn];
    auto *const t_nearest = new int[maxn];

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
            t_nearest[ncount++] = j;
          }
        }

        // if not nnn neighbors, order parameter = 0;

        if (ncount < nnn) {
          qn[0] = qn[1] = 0.0;
          continue;
        }

        // if nnn > 0, use only nearest nnn neighbors

        if (nnn > 0) {
          select2(nnn, ncount, t_distsq, t_nearest);
          ncount = nnn;
        }

        double usum = 0.0;
        double vsum = 0.0;

        for (int jj = 0; jj < ncount; jj++) {
          const int j = t_nearest[jj] & NEIGHMASK;
          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          double u, v;
          calc_qn_complex(delx, dely, u, v);
          usum += u;
          vsum += v;
        }
        qn[0] = usum / nnn;
        qn[1] = vsum / nnn;
      } else
        qn[0] = qn[1] = 0.0;
    }

    delete[] t_distsq;
    delete[] t_nearest;
  }
}

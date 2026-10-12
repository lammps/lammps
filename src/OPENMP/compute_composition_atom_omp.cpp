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

#include "compute_composition_atom_omp.h"

#include "atom.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCompositionAtomOMP::ComputeCompositionAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeCompositionAtom(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeCompositionAtom::compute_peratom().  The
   serial setup (array (re)allocation, neighbor-list build, zeroing) is
   unchanged; only the per-atom loop is threaded.  Each atom writes solely
   to its own output row (result[i]), so the result is bit-identical to
   the serial compute regardless of thread count or schedule.
------------------------------------------------------------------------- */

void ComputeCompositionAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow result array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(result);
    nmax = atom->nmax;
    memory->create(result, nmax, size_peratom_cols, "composition/atom:result");
    array_atom = result;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // zero the accumulators; atoms not in the group report only zeros

  memset(&result[0][0], 0, (size_t) nmax * size_peratom_cols * sizeof(double));

  // compute properties for each atom in group
  // use full neighbor list to count atoms less than cutoff

  const double *const *const x = atom->x;
  const int *const type = atom->type;
  const int *const mask = atom->mask;

  // get per-atom local compositions

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int ii = 0; ii < inum; ii++) {

    const int i = ilist[ii];

    if (mask[i] & groupbit) {

      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const int *const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      // i atom contribution

      int count = 1;

      const int itype = type[i];
      result[i][itype]++;

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;

        const int jtype = type[j];

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < cutsq) {
          count++;
          result[i][jtype]++;
        }
      }

      // total count of atoms found in sampled radius range

      result[i][0] = count;

      // local comp fractions per element

      const double lfac = 1.0 / count;
      for (int n = 1; n < size_peratom_cols; n++) result[i][n] *= lfac;
    }
  }
}

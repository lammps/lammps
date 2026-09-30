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

#include "compute_contact_atom_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeContactAtomOMP::ComputeContactAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeContactAtom(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeContactAtom::compute_peratom().  The serial
   setup (array (re)allocation, neighbor-list build) and the reverse
   communication are unchanged; the zeroing and tally loops are threaded.
   Since the tally scatters to BOTH atoms of a pair, the increments use
   atomic updates.  All tallied values are additions of exactly
   representable integers (+1.0), so the result is bit-identical to the
   serial compute regardless of thread count, schedule, or update order.
------------------------------------------------------------------------- */

void ComputeContactAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow contact array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(contact);
    nmax = atom->nmax;
    memory->create(contact, nmax, "contact/atom:contact");
    vector_atom = contact;
  }

  // invoke neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // compute number of contacts for each atom in group
  // contact if distance <= sum of radii
  // tally for both I and J

  const double *const *const x = atom->x;
  const double *const radius = atom->radius;
  const int *const mask = atom->mask;
  const int nall = atom->nlocal + atom->nghost;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < nall; i++) contact[i] = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];

    // only proceed if i is either part of the compute group or will contribute to contacts

    if (!(mask[i] & groupbit) && !(mask[i] & jgroupbit)) continue;

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const double radi = radius[i];
    const int *const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;

      // only tally for atoms in compute group (groupbit) if neighbor is in group2 (jgroupbit)

      const bool update_i_flag = (mask[i] & groupbit) && (mask[j] & jgroupbit);
      const bool update_j_flag = (mask[j] & groupbit) && (mask[i] & jgroupbit);
      if (!update_i_flag && !update_j_flag) continue;

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;
      const double radsum = radi + radius[j];
      const double radsumsq = radsum * radsum;
      if (rsq <= radsumsq) {
        if (update_i_flag) {
#if defined(_OPENMP)
#pragma omp atomic
#endif
          contact[i] += 1.0;
        }
        if (update_j_flag) {
#if defined(_OPENMP)
#pragma omp atomic
#endif
          contact[j] += 1.0;
        }
      }
    }
  }

  // communicate ghost atom counts between neighbor procs if necessary

  if (force->newton_pair) comm->reverse_comm(this);
}

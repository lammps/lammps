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

#include "pimd_partition.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "lammps.h"
#include "lmptype.h"
#include "universe.h"

#include <algorithm>
#include <vector>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   check per-atom type and group consistency across all partitions
------------------------------------------------------------------------- */

void PIMDUtils::check_atom_consistency(LAMMPS *lmp, const char *style, int groupbit)
{
  Atom *atom = lmp->atom;
  Comm *comm = lmp->comm;
  Universe *universe = lmp->universe;

  if (atom->natoms >= MAXSMALLINT)
    lmp->error->all(FLERR, "Fix {} atom type consistency check supports fewer than {} atoms", style,
                    MAXSMALLINT);
  const int natoms = static_cast<int>(atom->natoms);

  std::vector<tagint> tags(natoms);
  if (universe->iworld == 0) {
    std::vector<int> counts(comm->nprocs), displacements(comm->nprocs);
    MPI_Allgather(&atom->nlocal, 1, MPI_INT, counts.data(), 1, MPI_INT, lmp->world);
    for (int i = 1; i < comm->nprocs; i++) displacements[i] = displacements[i - 1] + counts[i - 1];
    MPI_Allgatherv(atom->tag, atom->nlocal, MPI_LMP_TAGINT, tags.data(), counts.data(),
                   displacements.data(), MPI_LMP_TAGINT, lmp->world);
    std::sort(tags.begin(), tags.end());
  }
  MPI_Bcast(tags.data(), natoms, MPI_LMP_TAGINT, 0, universe->uworld);

  std::vector<int> types(natoms, 0), members(natoms, 0);
  for (int i = 0; i < natoms; i++) {
    const int index = atom->map(tags[i]);
    if (index >= 0 && index < atom->nlocal) {
      types[i] = atom->type[index];
      members[i] = (atom->mask[index] & groupbit) ? 1 : 0;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, types.data(), natoms, MPI_INT, MPI_SUM, lmp->world);
  MPI_Allreduce(MPI_IN_PLACE, members.data(), natoms, MPI_INT, MPI_SUM, lmp->world);

  std::vector<int> reference_types(types);
  MPI_Bcast(reference_types.data(), natoms, MPI_INT, 0, universe->uworld);
  int mismatch = (types != reference_types) ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &mismatch, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (mismatch)
    lmp->error->all(
        FLERR, "Fix {} requires the same atom types for every atom ID in every partition", style);

  std::vector<int> reference_members(members);
  MPI_Bcast(reference_members.data(), natoms, MPI_INT, 0, universe->uworld);
  mismatch = (members != reference_members) ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &mismatch, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (mismatch)
    lmp->error->all(
        FLERR, "Fix {} requires the same group membership for every atom ID in every partition",
        style);
}

/* ----------------------------------------------------------------------
   collect requested owned-atom vectors within one bead world
------------------------------------------------------------------------- */

void PIMDUtils::collect_atom_vectors(LAMMPS *lmp, const char *style,
                                     const std::vector<tagint> &requested_tags, double **source,
                                     std::vector<double> &requested_values)
{
  Atom *atom = lmp->atom;
  Comm *comm = lmp->comm;
  Universe *universe = lmp->universe;

  const int me = comm->me;
  const int nprocs = comm->nprocs;
  const int nlocal = atom->nlocal;

  int too_many = (requested_tags.size() >= static_cast<std::size_t>(MAXSMALLINT)) ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &too_many, 1, MPI_INT, MPI_MAX, lmp->world);
  if (too_many)
    lmp->error->all(FLERR, "Fix {} atom vector collection supports fewer than {} requests", style,
                    MAXSMALLINT);

  const int nrequested = static_cast<int>(requested_tags.size());
  std::vector<int> counts(nprocs), displacements(nprocs);
  MPI_Allgather(&nrequested, 1, MPI_INT, counts.data(), 1, MPI_INT, lmp->world);
  std::size_t total_requested_size = 0;
  for (int rank = 0; rank < nprocs; rank++) {
    displacements[rank] = static_cast<int>(total_requested_size);
    total_requested_size += counts[rank];
    if (total_requested_size >= static_cast<std::size_t>(MAXSMALLINT))
      lmp->error->all(FLERR, "Fix {} atom vector collection supports fewer than {} requests", style,
                      MAXSMALLINT);
  }
  const int total_requested = static_cast<int>(total_requested_size);

  std::vector<tagint> all_tags(total_requested);
  MPI_Allgatherv(requested_tags.data(), nrequested, MPI_LMP_TAGINT, all_tags.data(), counts.data(),
                 displacements.data(), MPI_LMP_TAGINT, lmp->world);

  std::vector<int> owners(total_requested, 0);
  std::vector<double> all_values(3 * static_cast<std::size_t>(total_requested), 0.0);
  for (int i = 0; i < total_requested; i++) {
    const int index = atom->map(all_tags[i]);
    if (index >= 0 && index < nlocal) {
      owners[i] = 1;
      all_values[3 * i] = source[index][0];
      all_values[3 * i + 1] = source[index][1];
      all_values[3 * i + 2] = source[index][2];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, owners.data(), total_requested, MPI_INT, MPI_SUM, lmp->world);
  for (std::size_t offset = 0; offset < all_values.size();) {
    const int count = static_cast<int>(
        std::min(all_values.size() - offset, static_cast<std::size_t>(MAXSMALLINT)));
    MPI_Allreduce(MPI_IN_PLACE, all_values.data() + offset, count, MPI_DOUBLE, MPI_SUM, lmp->world);
    offset += count;
  }

  requested_values.resize(3 * requested_tags.size());
  const int first = displacements[me];
  for (int i = 0; i < nrequested; i++) {
    const int index = first + i;
    if (owners[index] != 1)
      lmp->error->universe_one(FLERR,
                               fmt::format("Fix {} found {} owners for atom ID {} in partition {}",
                                           style, owners[index], requested_tags[i],
                                           universe->iworld));
    requested_values[3 * i] = all_values[3 * index];
    requested_values[3 * i + 1] = all_values[3 * index + 1];
    requested_values[3 * i + 2] = all_values[3 * index + 2];
  }
}

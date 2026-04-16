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

/* ----------------------------------------------------------------------
   fix msevb transfer:
   Potential energy gathering and permanent topology transfer logic.
------------------------------------------------------------------------- */

#include "fix_msevb.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "memory.h"
#include "neighbor.h"
#include "special.h"
#include "universe.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   Broadcast all ref_* topology arrays from partition 0 to all others.
   Requires snapshot_reference_topology() to have been called first.
---------------------------------------------------------------------- */

void FixMSEVB::broadcast_reference_topology() {
  const tagint sz = ref_maxtag + 1;
  MPI_Bcast(ref_type, sz, MPI_INT, 0, samerank);
  MPI_Bcast(ref_charge, sz, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(ref_molecule, sz, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_num_bond, sz, MPI_INT, 0, samerank);
  MPI_Bcast(ref_bond_type_flat, sz * ref_bond_per_atom, MPI_INT, 0, samerank);
  MPI_Bcast(ref_bond_atom_flat, sz * ref_bond_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_num_angle, sz, MPI_INT, 0, samerank);
  MPI_Bcast(ref_angle_type_flat, sz * ref_angle_per_atom, MPI_INT, 0, samerank);
  MPI_Bcast(ref_angle_atom1_flat, sz * ref_angle_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_angle_atom2_flat, sz * ref_angle_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_angle_atom3_flat, sz * ref_angle_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_num_dihedral, sz, MPI_INT, 0, samerank);
  MPI_Bcast(ref_dihedral_type_flat, sz * ref_dihedral_per_atom, MPI_INT, 0, samerank);
  MPI_Bcast(ref_dihedral_atom1_flat, sz * ref_dihedral_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_dihedral_atom2_flat, sz * ref_dihedral_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_dihedral_atom3_flat, sz * ref_dihedral_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_dihedral_atom4_flat, sz * ref_dihedral_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_num_improper, sz, MPI_INT, 0, samerank);
  MPI_Bcast(ref_improper_type_flat, sz * ref_improper_per_atom, MPI_INT, 0, samerank);
  MPI_Bcast(ref_improper_atom1_flat, sz * ref_improper_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_improper_atom2_flat, sz * ref_improper_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_improper_atom3_flat, sz * ref_improper_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_improper_atom4_flat, sz * ref_improper_per_atom, MPI_LMP_TAGINT, 0, samerank);
  MPI_Bcast(ref_nspecial_flat, sz * 3, MPI_INT, 0, samerank);
  MPI_Bcast(ref_special_flat, sz * ref_maxspecial, MPI_LMP_TAGINT, 0, samerank);
}

/* ----------------------------------------------------------------------
   Gather potential energy from each partition via MPI_Allreduce.
---------------------------------------------------------------------- */

void FixMSEVB::gather_potential_energies() {
  const int np = npartitions;
  const double scalefac = 1.0 / comm->nprocs;

  double t0 = MPI_Wtime();
  for (int i = 0; i < np; i++)
    my_epot[i] = 0.0;
  my_epot[ipartition] = scalefac * pe->compute_scalar();
  double t1 = MPI_Wtime();
  prof_pe_compute += t1 - t0;

  // Start non-blocking allreduce — completion waited in post_force()
  // after gather_forces() overlaps with the MPI progress.
  t0 = MPI_Wtime();
  MPI_Iallreduce(my_epot, epot, npartitions, MPI_DOUBLE, MPI_SUM,
                 universe->uworld, &epot_request);
  t1 = MPI_Wtime();
  prof_pe_allreduce += t1 - t0;

  pe->addstep(update->ntimestep + 1);
}

/* ----------------------------------------------------------------------
   Permanent transfer: if the most probable EVB state is not the
   reference (state 0), commit the topology change to the reference.
   Matches Python MuStaRD behavior — stateless, no threshold or cooldown.
---------------------------------------------------------------------- */

bool FixMSEVB::do_permanent_transfer(int &out_max_state, double &out_max_amp) {
  if (!reaction_enabled || nsites == 0) {
    out_max_state = 0;
    out_max_amp   = (amplitudes ? amplitudes[0] : 0.0);
    return false;
  }

  // Find state with maximum amplitude
  int max_state = 0;
  double max_amp = amplitudes[0];
  for (int k = 1; k < nstates; k++) {
    if (amplitudes[k] > max_amp) {
      max_amp = amplitudes[k];
      max_state = k;
    }
  }

  if (max_state == 0) {
    out_max_state = 0;
    out_max_amp   = max_amp;
    return false;
  }
  int sk = max_state - 1;

  // Commit: apply every step in the chain to partition 0's topology.
  if (ipartition == 0) {
    if (sites[sk].n_components > 0) {
      // Product state: apply every depth of every component's chain.
      for (int ci = 0; ci < sites[sk].n_components; ci++) {
        int comp = sites[sk].components[ci];
        const int comp_len = sites[comp].chain_len;
        for (int d = 0; d < comp_len; d++) {
          const ReactionDef &rxn_c =
              rxndefs[chain_rxn_flat[comp * max_shells + d]];
          const tagint *glove_c =
              chain_glove_flat + (comp * max_shells + d) * glove_nmax;
          tagint tX = chain_X_flat[comp * max_shells + d];
          tagint tH = chain_H_flat[comp * max_shells + d];
          tagint tY = chain_Y_flat[comp * max_shells + d];
          apply_state_change(glove_c, atom->map(tX), atom->map(tH),
                             atom->map(tY), tX, tH, tY, rxn_c);
          if (d < comp_len - 1)
            comm->forward_comm(this);
        }
        if (ci < sites[sk].n_components - 1)
          comm->forward_comm(this);
      }
    } else {
      const int chain_len = sites[sk].chain_len;

      for (int d = 0; d < chain_len; d++) {
        const ReactionDef &rxn_d =
            rxndefs[chain_rxn_flat[sk * max_shells + d]];
        const tagint *glove_d =
            chain_glove_flat + (sk * max_shells + d) * glove_nmax;
        tagint tX = chain_X_flat[sk * max_shells + d];
        tagint tH = chain_H_flat[sk * max_shells + d];
        tagint tY = chain_Y_flat[sk * max_shells + d];
        apply_state_change(glove_d, atom->map(tX), atom->map(tH),
                           atom->map(tY), tX, tH, tY, rxn_d);
      }
    }
  }

  out_max_state = max_state;
  out_max_amp   = max_amp;

  if (universe->me == 0) {
    auto msg = fmt::format("MSEVB: permanent transfer at step {} "
                           "(site {}, amplitude {:.4f})\n",
                           update->ntimestep, sk, max_amp);
    utils::logmesg(lmp, msg);
  }

  // Partition 0 now has the new reference topology (apply_state_change
  // above). Non-0 partitions still have their per-state modified topology.
  // Snapshot partition 0's new topology and broadcast to all.
  snapshot_reference_topology();

  broadcast_reference_topology();

  // Restore the new reference topology on ALL partitions
  restore_reference_topology();

  // Rebuild special bonds from the new bond topology, then do a full
  // reneighbor sequence so that atoms are in the correct subdomains
  // and ghost atoms have the new topology. This is essential for
  // multi-rank partitions where atom migration may be needed.
  {
    Special special(lmp);
    special.build(true);
    // Notify the Kokkos DualView that special bonds were updated on the host,
    // and push them to device.  Without this, sync_before_neighbor_build()
    // sees no host modification and skips the host→device push, so
    // neighbor->build() bakes stale pre-transfer specials into the GPU
    // neighbor list.  (CPU no-op; critical for msevb/kk.)
    modified_topology_on_host();
    if (domain->triclinic)
      domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    comm->borders();
    if (domain->triclinic)
      domain->lamda2x(atom->nlocal + atom->nghost);
    sync_before_neighbor_build();
    neighbor->build(1);
  }

  // Re-snapshot so the reference includes updated specials
  snapshot_reference_topology();

  // Sync velocities from partition 0 to other partitions.
  // samerank connects same-rank procs across partitions (which have
  // the same nlocal and atom ordering), so this is safe.
  if (atom->nlocal > 0)
    MPI_Bcast(&atom->v[0][0], 3 * atom->nlocal, MPI_DOUBLE, 0, samerank);
  return true;
}

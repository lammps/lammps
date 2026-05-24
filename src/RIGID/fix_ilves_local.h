/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(ilves/local,FixIlvesLocal);
// clang-format on
#else

#ifndef LMP_FIX_ILVES_LOCAL_H
#define LMP_FIX_ILVES_LOCAL_H

#include "fix_ilves.h"

namespace LAMMPS_NS {

// Local-topology ILVES variant.  Builds constraint clusters purely from
// each rank's local bond / angle storage, with no global topology gather.
// Per-rank memory therefore scales with the number of local atoms rather
// than with the total system size.
//
// Restrictions:
//   * `newton off bond` is required (each bond must be stored on both
//     endpoints).
//   * Every constraint cluster must fit within a single subdomain plus
//     the communication cutoff -- i.e. all atoms of a cluster must be
//     either local or available as ghosts on the rank that owns the
//     cluster's center.  If validation fails, the fix errors out at init
//     with instructions to switch to fix ilves/global or extend
//     `comm_modify cutoff`.
//   * Cluster topology must be "star-shaped" (one central atom with k
//     bonds to leaves, k >= 1) or a 1+1 pair (single bond).  Multi-hop
//     topologies such as a 3+-bond polymer backbone are rejected; use
//     fix ilves/global for those.

class FixIlvesLocal : public FixIlves {
 public:
  FixIlvesLocal(class LAMMPS *, int, char **);
  ~FixIlvesLocal() override;

  // override per-atom array lifecycle so the reverse-comm buffer grows
  // alongside the base class's per-atom arrays
  double memory_usage() override;
  void grow_arrays(int) override;

  // reverse_comm sends per-atom constraint contributions (forces in
  // post_force, velocity deltas in end_of_step) from ghost atoms on the
  // cluster owner back to the corresponding local atoms on the owner of
  // each ghost.  Both endpoints pack/unpack the same 3-double payload.
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  void init_topology() override;
  void build_constraint_list() override;
  bool bond_is_constrained(tagint ta, tagint tb) override;

  // Cluster-owner force / velocity application path with reverse_comm.
  // Overrides the base-class redundant-solve implementations.
  void apply_constraint_forces(int vflag) override;
  void correct_coordinates(int vflag) override;
  void correct_velocities() override;

  // No-op: the local variant maintains ghost xshake locally and does not
  // forward_comm during Newton iteration.
  void sync_xshake() override {}

  // The cluster owner updates ghost atoms' xshake locally because the
  // ghost-owning ranks don't process this cluster.
  bool update_ghost_xshake() const override { return true; }

 private:
  // Per-atom buffer used to accumulate this fix's contribution at ghost
  // atoms during one of {force, position, velocity} application phases,
  // then reverse_comm'd so the owning rank receives the accumulated
  // contribution at its local index.  Reused across phases because the
  // three never overlap in time.
  double **rbuf;
  int maxbuf;

  // count the number of locally-stored constrained bonds anchored at
  // local atom i (i must be < nlocal).
  int count_constrained_bonds_local(int i);

  // build angle_distance[] from local angle storage with an MPI_Allreduce
  // to share flanking-bond-type information across ranks that may not
  // each see every angle type locally.
  void compute_angle_distances_local();

  void grow_rbuf(int);
};

}    // namespace LAMMPS_NS

#endif
#endif

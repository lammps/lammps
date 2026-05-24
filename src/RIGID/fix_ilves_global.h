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
FixStyle(ilves/global,FixIlvesGlobal);
// clang-format on
#else

#ifndef LMP_FIX_ILVES_GLOBAL_H
#define LMP_FIX_ILVES_GLOBAL_H

#include "fix_ilves.h"

#include <unordered_map>
#include <vector>

namespace LAMMPS_NS {

// Global-topology ILVES variant.  Gathers the complete bond/angle table
// onto every MPI rank once at init via MPI_Allgatherv, then builds the
// per-rank constraint list to include every cluster that intersects at
// least one local atom -- even constraints between two ghost atoms.
// Suitable for constraint clusters that span multiple subdomains (e.g.
// an all-backbone-constrained polymer); the cost is that the replicated
// topology grows with system size rather than with subdomain size, so
// per-rank memory does not shrink with rank count.

class FixIlvesGlobal : public FixIlves {
 public:
  FixIlvesGlobal(class LAMMPS *, int, char **);
  ~FixIlvesGlobal() override = default;

 protected:
  void init_topology() override;
  void build_constraint_list() override;
  bool bond_is_constrained(tagint ta, tagint tb) override;

 private:
  // Replicated bond / angle tables, gathered once at init.
  // Bonds: (lower-tag, higher-tag, type), sorted by (a, b) and deduped.
  // Angles: (atom1, mid, atom3, type), sorted by middle then outer; only
  // angle types listed in angle_flag[] are gathered.
  std::vector<tagint> gb_a, gb_b;
  std::vector<int> gb_type;
  std::vector<tagint> ga1, ga2, ga3;
  std::vector<int> ga_type;
  bool global_topology_ready;

  // global cluster-ID for every tag involved in any bond/angle.
  // Maps a tag to its cluster's representative tag.  Sparse: only tags
  // that participate in at least one constrained bond or angle appear in
  // the map.  This drops the per-rank memory from ~8*natoms (the dense
  // std::vector<tagint> previously used) to ~24-40 bytes per involved
  // tag (typical std::unordered_map overhead), which is the only viable
  // shape for partial-constraint systems where the involved subset is
  // much smaller than natoms (e.g. waters in a large solvent box).
  std::unordered_map<tagint, tagint> tag_cluster;

  void gather_global_topology();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* -*- c++ -*- ----------------------------------------------------------
   MSEVB superimpose: BFS graph-isomorphism for template matching.
   Derived from REACTER's superimpose algorithm (Gissinger et al.).
------------------------------------------------------------------------- */

#ifndef LMP_FIX_MSEVB_SUPERIMPOSE_H
#define LMP_FIX_MSEVB_SUPERIMPOSE_H

#include "lammps.h"
#include "molecule.h"

#include <unordered_map>
#include <vector>

namespace LAMMPS_NS {

// VirtualTopo: lightweight overlay for multi-shell detection.
// When populated, these maps override atom->type[] and atom->special[]
// for atoms that appear as keys.  Keys are global atom tags.
struct VirtualTopo {
  std::unordered_map<tagint, int>                type_override;
  std::unordered_map<tagint, std::vector<tagint>> bond_override;
};

// RefTopo: reference topology snapshot used as fallback for ghost/off-rank atoms.
// Populated from FixMSEVB's ref_* arrays after snapshot_reference_topology().
// Allows superimpose to traverse atoms that are not local to the calling rank.
//
// Uses the special list (symmetric 1-2 neighbors) rather than the one-sided
// bond_atom list, because LAMMPS stores each bond on only one atom but the
// special list is symmetric (both bonded atoms list each other).
struct RefTopo {
  const int    *type;           // [0..maxtag]: atom type by global tag
  const int    *nspecial_flat;  // [tag*3 + 0]: count of 1-2 bonded neighbors
  const tagint *special_flat;   // [tag*maxspecial + k]: k-th 1-2 neighbor tag
  int           maxspecial;     // stride for special_flat
  tagint        maxtag;         // maximum valid tag index
};

// msevb_superimpose: attempt to match molecule template mol onto the real
// system starting from the seed pair (tag_H → ibonding, tag_Y → jbonding).
//
// Parameters:
//   lmp       — LAMMPS instance (provides atom->type, atom->special, etc.)
//   mol       — pre-reaction molecule template
//   is_edge   — [mol->natoms]: 1 if atom is an edge/boundary atom (not matched)
//   ibonding  — 0-based index in mol of the H (initiating) atom
//   jbonding  — 0-based index in mol of the Y (accepting) atom
//   tag_H     — global tag of the real H atom
//   tag_Y     — global tag of the real Y atom
//   glove     — output: [mol->natoms] global tags of matched real atoms (0=unmatched)
//   vtopo     — optional virtual topology overlay for multi-shell detection
//   ref       — optional reference topology for ghost/off-rank atom fallback
//
// Returns true on success (all non-edge atoms assigned uniquely), false otherwise.
bool msevb_superimpose(LAMMPS *lmp, Molecule *mol, const int *is_edge,
                       int ibonding, int jbonding,
                       tagint tag_H, tagint tag_Y, tagint *glove,
                       const VirtualTopo *vtopo = nullptr,
                       const RefTopo     *ref   = nullptr);

} // namespace LAMMPS_NS

#endif

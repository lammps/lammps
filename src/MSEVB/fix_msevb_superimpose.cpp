/* ----------------------------------------------------------------------
   MS-EVB superimpose: stripped-down graph isomorphism for template matching.
   Derived from REACTER's superimpose algorithm (Gissinger et al.).
------------------------------------------------------------------------- */

#include "fix_msevb_superimpose.h"

#include "atom.h"
#include "molecule.h"

#include <cstring>
#include <vector>
#include <unordered_set>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   msevb_superimpose

   Molecule arrays use 0-based atom indexing (0 .. natoms-1):
     mol->type[i]            — atom type of template atom i
     mol->nspecial[i][0]     — number of 1-2 bonded neighbors of atom i
     mol->special[i][j]      — 1-based atom ID of the j-th 1-2 neighbor;
                               subtract 1 to get 0-based glove index.

   glove[i] holds the real global tag corresponding to template atom i.
   A value of 0 means "not yet assigned".

   vtopo (optional): if non-null, its type_override and bond_override maps
   take precedence over atom->type[] and atom->special[] respectively.
   bond_override values are vectors of global tags (not local indices).

   ref (optional): reference topology snapshot used as fallback when an atom
   is not local (ghost or off-rank).  Allows superimpose to work in multi-rank
   partitions where some template atoms are owned by other ranks.
---------------------------------------------------------------------- */

bool LAMMPS_NS::msevb_superimpose(LAMMPS *lmp, Molecule *mol,
                                   const int *is_edge, int ibonding,
                                   int jbonding, tagint tag_H, tagint tag_Y,
                                   tagint *glove,
                                   const VirtualTopo *vtopo,
                                   const RefTopo     *ref) {
  const int natoms = mol->natoms;
  memset(glove, 0, natoms * sizeof(tagint));
  glove[ibonding] = tag_H;
  glove[jbonding] = tag_Y;

  std::unordered_set<tagint> assigned;
  assigned.insert(tag_H);
  assigned.insert(tag_Y);

  Atom *atom = lmp->atom;
  int *atype = atom->type;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  const int nlocal = atom->nlocal;

  // --- Helpers -----------------------------------------------------------

  // Number of 1-2 template neighbors of template atom i (0-based).
  auto tmpl_nn = [&](int i) { return mol->nspecial[i][0]; };
  // j-th 1-2 template neighbor of template atom i, as 0-based index.
  auto tmpl_neigh = [&](int i, int j) {
    return (int)(mol->special[i][j] - 1); // convert from 1-based
  };

  // Check whether template atom ti has any unassigned non-edge neighbor.
  auto has_unassigned_neigh = [&](int ti) -> bool {
    for (int j = 0; j < tmpl_nn(ti); j++) {
      int k = tmpl_neigh(ti, j);
      if (!is_edge[k] && glove[k] == 0)
        return true;
    }
    return false;
  };

  // Get type of real atom with global tag r.
  // Priority: vtopo override > local atom->type > ref snapshot.
  auto real_type = [&](tagint r) -> int {
    if (vtopo) {
      auto it = vtopo->type_override.find(r);
      if (it != vtopo->type_override.end())
        return it->second;
    }
    int local_r = atom->map(r);
    if (local_r >= 0 && local_r < nlocal)
      return atype[local_r];
    // Fall back to ref snapshot for ghost/off-rank atoms.
    if (ref && r >= 0 && r <= ref->maxtag)
      return ref->type[r];
    return 0;
  };

  // Get the bond-partner tag list for real atom with global tag r.
  // Priority: vtopo override > local atom->special > ref snapshot.
  // Returns an empty vector if no data is available.
  auto get_bond_tags = [&](tagint r) -> std::vector<tagint> {
    // vtopo override takes priority.
    if (vtopo) {
      auto vit = vtopo->bond_override.find(r);
      if (vit != vtopo->bond_override.end())
        return vit->second;
    }
    // Local atom: use LAMMPS special list.
    int local_r = atom->map(r);
    if (local_r >= 0 && local_r < nlocal) {
      int ns = nspecial[local_r][0];
      std::vector<tagint> tags(ns);
      for (int s = 0; s < ns; s++)
        tags[s] = special[local_r][s];
      return tags;
    }
    // Ghost / off-rank: use reference special list (symmetric 1-2 neighbors).
    if (ref && r >= 0 && r <= ref->maxtag) {
      int nb = ref->nspecial_flat[r * 3 + 0]; // count of 1-2 bonded neighbors
      std::vector<tagint> tags(nb);
      for (int b = 0; b < nb; b++)
        tags[b] = ref->special_flat[(tagint)r * ref->maxspecial + b];
      return tags;
    }
    return {};
  };

  // True if assigning glove[tneigh] = r is consistent with the current glove:
  // every already-assigned template neighbor k of tneigh must appear in
  // the bond-partner list of r.
  auto compatible = [&](int tneigh, tagint r) -> bool {
    std::vector<tagint> nbrs = get_bond_tags(r);
    if (nbrs.empty())
      return false; // can't verify; fail safely

    for (int j = 0; j < tmpl_nn(tneigh); j++) {
      int k = tmpl_neigh(tneigh, j);
      if (glove[k] == 0)
        continue; // not yet assigned; skip
      bool found = false;
      for (tagint nb_tag : nbrs) {
        if (nb_tag == glove[k]) {
          found = true;
          break;
        }
      }
      if (!found)
        return false;
    }
    return true;
  };

  // --- BFS queue of "pioneer" template atoms ----------------------------
  // A pioneer is assigned and has at least one unassigned non-edge neighbor.

  std::vector<int> queue;
  queue.reserve(natoms);

  if (has_unassigned_neigh(ibonding))
    queue.push_back(ibonding);
  if (has_unassigned_neigh(jbonding))
    queue.push_back(jbonding);

  // Process pioneers.
  for (size_t qi = 0; qi < queue.size(); qi++) {
    const int pion = queue[qi];
    const tagint pion_tag = glove[pion];

    // Get bond-partner tags for this pioneer's real atom.
    std::vector<tagint> pion_nbrs = get_bond_tags(pion_tag);
    if (pion_nbrs.empty())
      return false; // no connectivity data available for this pioneer

    // For each unassigned non-edge template neighbor of this pioneer:
    for (int j = 0; j < tmpl_nn(pion); j++) {
      const int tneigh = tmpl_neigh(pion, j);
      if (is_edge[tneigh])
        continue;
      if (glove[tneigh] != 0)
        continue; // already assigned

      const int ttype = mol->type[tneigh];

      // Find a valid real atom among the pioneer's bonded neighbors that:
      //   (a) has the right type (with vtopo/ref override if applicable),
      //   (b) is not already in the glove, and
      //   (c) is compatible with all other already-assigned template neighbors.
      // Take the first valid candidate (handles symmetric atoms).
      tagint candidate = 0;

      for (tagint r : pion_nbrs) {
        int rt = real_type(r);
        if (rt != ttype)
          continue;
        if (assigned.count(r))
          continue;
        if (!compatible(tneigh, r))
          continue;

        candidate = r;
        break; // take first valid symmetric candidate
      }

      if (candidate == 0)
        return false; // no match found

      glove[tneigh] = candidate;
      assigned.insert(candidate);

      // If tneigh has unassigned non-edge neighbors, make it a pioneer.
      if (has_unassigned_neigh(tneigh)) {
        bool already_queued = false;
        for (size_t q2 = 0; q2 < queue.size(); q2++) {
          if (queue[q2] == tneigh) {
            already_queued = true;
            break;
          }
        }
        if (!already_queued)
          queue.push_back(tneigh);
      }
    }
  }

  // Verify every non-edge template atom is assigned.
  for (int i = 0; i < natoms; i++)
    if (!is_edge[i] && glove[i] == 0)
      return false;

  return true;
}

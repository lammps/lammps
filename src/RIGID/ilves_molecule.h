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

/* ----------------------------------------------------------------------
   ILVES constraint solver: constraint topology (atom graph + bond graph).
   Ported from GROMACS 2021 ILVES (LGPL-2.1), src/gromacs/mdlib/molecule.{cpp,h}.
   The constructor is adapted to take a LAMMPS-built constraint list (atom
   index pairs + target distances) instead of the GROMACS InteractionDefinitions
   / domain-decomposition structures; the graph construction is unchanged.
   See ilves_graph.h for full attribution.
------------------------------------------------------------------------- */

#ifndef LMP_ILVES_MOLECULE_H
#define LMP_ILVES_MOLECULE_H

#include "ilves_graph.h"

#include <vector>

/**
 * The Molecule structure contains all the information about the bonds and atoms
 * of the molecule required by ILVES constraint solver. A Molecule structure can
 * contain several molecules, with out any connection between them.
 *
 */

namespace LAMMPS_NS::ILVES {

class Molecule {
 public:
  struct Atoms {
    // invmass[i] is 1/mass of the ith atom (borrowed; owned by the caller).
    // The atom graph and atom count are construction-only and not retained.
    const double *invmass;
  } atoms;

  struct Bonds {
    int num;        // The number of bonds in the molecule.
    Graph graph;    // The bond graph of the molecule.

    // GEOMETRY indices of the two atoms of each bond: the nearest periodic
    // image (local or ghost) whose coordinates give the correct bond vector
    // by raw subtraction at any box size.  Used for positions, increments,
    // and the virial.  Example: bond 0 joins atom1[0] and atom2[0].
    std::vector<int> atom1;
    std::vector<int> atom2;

    // NODE (canonical owner) indices of the same two atoms: an atom and its
    // periodic ghost image share one node id, so a bond that wraps a
    // periodic boundary stays a single edge and the constraint graph does
    // not fragment.  Used to build the graph and to detect the shared atom
    // when assembling the Jacobian (see Ilves::make_weights).  node1[k]
    // pairs with atom1[k] (same physical atom), node2[k] with atom2[k].
    std::vector<int> node1;
    std::vector<int> node2;

    // The bond length squared of each bond.
    std::vector<double> sigma2;
    std::vector<double> invsigma2;
  } bonds;

  /**
     * Default constructor is deleted to avoid errors.
     */
  Molecule() = delete;

  /**
     * Constructs a Molecule that can be used by the ILVES constraint solver,
     * from a LAMMPS-built constraint list.
     *
     * The number of constraints is catom1.size().
     * @param catom1 catom1[k]/catom2[k] are the GEOMETRY indices (nearest
     * periodic image, into the position / force / invmass arrays) of the two
     * atoms joined by constraint k.
     * @param catom2 See catom1.
     * @param cnode1 cnode1[k]/cnode2[k] are the NODE (canonical owner) indices
     * of the same two atoms, used only for graph connectivity and shared-atom
     * detection; cnode1[k] is the owner of catom1[k], cnode2[k] of catom2[k].
     * @param cnode2 See cnode1.
     * @param cdist cdist[k] is the target length of constraint k.
     * @param invmass The inverse mass of each atom.  The Molecule keeps a
     * pointer into this vector, so it must be kept alive (and not reallocated)
     * for the Molecule's lifetime.
     */
  Molecule(const std::vector<int> &catom1, const std::vector<int> &catom2,
           const std::vector<int> &cnode1, const std::vector<int> &cnode2,
           const std::vector<double> &cdist, const std::vector<double> &invmass);

  /**
     * Renumber the data of the Bonds structure given a permutation.
     * The permutation is given as in MATLAB. Example:
     *  p = [2, 1, 0] Means that
     *  Old position 2 is now position 0
     *  Old position 1 is now position 1
     *  Old position 0 is now position 2
     *
     * @param perm The permutation in MATLAB format.
     */
  void renumber_bonds(const std::vector<int> &perm);

  /**
     * Estimate the memory used by the constraint topology (atom and bond
     * graphs plus the per-bond index and distance arrays).  The borrowed
     * inverse-mass array is owned by the caller and is not counted.
     *
     * @return The size of the topology storage in bytes.
     */
  [[nodiscard]] double memory_usage() const;
};

}    // namespace LAMMPS_NS::ILVES

#endif

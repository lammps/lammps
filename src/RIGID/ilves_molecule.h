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

namespace LAMMPS_NS {
namespace ILVES {

class Molecule {
public:
    struct Atoms {
        int num;       // The number of atoms in the molecule. Only atoms with
                       // constraints.
        Graph graph;   // the atomic graph of the molecule. Only atoms with
                       // constraints.

        const double *invmass;   // invmass[i] is 1/mass of the ith atom
    } atoms;

    struct Bonds {
        int num;       // The number of bonds in the molecule.
        Graph graph;   // The bond graph of the molecule.

        // GROMACS atom indices of the two atoms of each bond. Example:
        // Bond 0 connects atom1[0] and atom2[0].
        std::vector<int> atom1;
        std::vector<int> atom2;

        // Local atom indices of the two atoms of each bond.
        // The local atom indices correspond to the numbering of
        // atoms_data.graph. Example: Bond 0 connects latom1[0] and latom2[0].
        std::vector<int> latom1;
        std::vector<int> latom2;

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
     * @param nbonds Number of constraints (bonds).
     * @param catom1 catom1[k]/catom2[k] are the indices (into the position /
     * force / invmass arrays) of the two atoms joined by constraint k.
     * @param catom2 See catom1.
     * @param cdist cdist[k] is the target length of constraint k.
     * @param invmass A pointer to the array of inverse mass of each atom.
     * The Molecule keeps a pointer to this array, so it must be kept alive.
     */
    Molecule(int nbonds, const int *catom1, const int *catom2, const double *cdist,
             const double *invmass);

    /**
     * Renumber the data of the Bonds structure given a permutation.
     * The permutation is given as in MATLAB. Example:
     *  p = [2, 1, 0] Means that
     *  Old position 2 is now position 0
     *  Old position 1 is now position 1
     *  Old position 0 is now position 2
     *
     * @param perm The permutation in MATLAB format.
     * @param renumber_graph If true, the bond graph is also renumbered. The
     * bond graph is not renumbered otherwise.
     */
    void renumber_bonds(const std::vector<int> &perm, bool renumber_graph);

    /**
     * Estimate the memory used by the constraint topology (atom and bond
     * graphs plus the per-bond index and distance arrays).  The borrowed
     * inverse-mass array is owned by the caller and is not counted.
     *
     * @return The size of the topology storage in bytes.
     */
    double memory_usage() const;
};

}    // namespace ILVES
}    // namespace LAMMPS_NS

#endif

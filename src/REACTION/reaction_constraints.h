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
   Contributing Author: Jacob Gissinger (jgissing@stevens.edu)
------------------------------------------------------------------------- */

#ifndef LMP_REACTION_CONSTRAINTS_H
#define LMP_REACTION_CONSTRAINTS_H

#include "pointers.h"    // IWYU pragma: export
#include "reaction.h"

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace LAMMPS_NS {

class ReactionConstraints : protected Pointers {
public:
  ReactionConstraints(LAMMPS *lmp);
  ~ReactionConstraints();

  int check(Reaction &, std::vector<tagint> &);
  void customvarnames(std::vector<Reaction> &);                     // get per-atom variables names used by custom constraint
  double get_temperature(std::vector<tagint> &);
  void get_customvars(int);                                // evaluate local values for variables names used by custom constraint

  int ncustomvars;
  int atoms2bondflag;                                      // 1 if atoms2bond map has been populated on this timestep
  double **vvec;

private:
  int nrxnfunction;
  std::vector<std::string> rxnfunclist;                    // lists current special rxn function
  std::vector<int> peratomflag;                            // 1 if special rxn function uses per-atom variable (vs. per-bond)

  std::vector<std::string> customvarstrs;
  int nvvec;                                               // per-atom vector to store custom constraint atom-style variable values

  std::map<std::set<tagint>, int> atoms2bond;              // maps atom pair to index of local bond array
  class Compute *cperbond;                                 // pointer to 'compute bond/local' used by custom constraint ('rxnbond' function)

  void get_IDcoords(Reaction::Constraint::IDType, int, double *, Molecule *, std::vector<tagint> &);
  bool custom_constraint(const std::string &, Reaction &, std::vector<tagint> &);
  double rxnfunction(const std::string &, const std::string &, const std::string &, Molecule *, std::vector<tagint> &);
  void get_atoms2bond(int);
};

}    // namespace LAMMPS_NS

#endif

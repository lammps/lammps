/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(delete_atoms,DeleteAtoms);
// clang-format on
#else

#ifndef LMP_DELETE_ATOMS_H
#define LMP_DELETE_ATOMS_H

#include "command.h"
#include <map>

namespace LAMMPS_NS {

class DeleteAtoms : public Command {
 public:
  DeleteAtoms(class LAMMPS *);
  void command(int, char **) override;

 private:
  int *dlist;
  int allflag, compress_flag, bond_flag, mol_flag;
  std::map<tagint, int> *hash;

  void delete_group(int, char **);
  void delete_region(int, char **);
  void delete_overlap(int, char **);
  void delete_random(int, char **);
  void delete_variable(int, char **);

  void delete_bond();
  void delete_molecule();
  void recount_topology();
  void options(int, char **);

  inline int sbmask(int j) const { return j >> SBBITS & 3; }

  // callback functions for ring communication

  static void bondring(int, char *, void *);
  static void molring(int, char *, void *);
};

}    // namespace LAMMPS_NS

#endif
#endif

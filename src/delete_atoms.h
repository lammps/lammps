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
  void command(int, char **);

 private:
  int *dlist;
  int allflag, compress_flag, bond_flag, mol_flag;
  std::map<tagint, int> *hash;

  void delete_group(int, char **);
  void delete_region(int, char **);
  void delete_overlap(int, char **);
  void delete_porosity(int, char **);

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

/* ERROR/WARNING messages:

E: Delete_atoms command before simulation box is defined

The delete_atoms command cannot be used before a read_data,
read_restart, or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use delete_atoms unless atoms have IDs

Your atoms do not have IDs, so the delete_atoms command cannot be
used.

W: Attempting to delete atoms in rigid bodies

UNDOCUMENTED

W: Ignoring 'compress yes' for molecular system

UNDOCUMENTED

E: Could not find delete_atoms group ID

Group ID used in the delete_atoms command does not exist.

E: Could not find delete_atoms region ID

Region ID used in the delete_atoms command does not exist.

E: Delete_atoms requires a pair style be defined

This is because atom deletion within a cutoff uses a pairwise
neighbor list.

E: Delete_atoms cutoff > max neighbor cutoff

Can only delete atoms in atom pairs that will be in neighbor list.

W: Delete_atoms cutoff > minimum neighbor cutoff

This means atom pairs for some atom types may not be in the neighbor
list and thus an atom in that pair cannot be deleted.

E: Cannot delete_atoms bond yes for non-molecular systems

Self-explanatory.

E: Cannot use delete_atoms bond yes with atom_style template

This is because the bonds for that atom style are hardwired in the
molecule template.

E: Delete_atoms mol yes requires atom attribute molecule

Cannot use this option with a non-molecular system.

*/

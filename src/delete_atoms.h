/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(delete_atoms,DeleteAtoms)

#else

#ifndef LMP_DELETE_ATOMS_H
#define LMP_DELETE_ATOMS_H

#include "pointers.h"
#include <map>

namespace LAMMPS_NS {

class DeleteAtoms : protected Pointers {
 public:
  DeleteAtoms(class LAMMPS *);
  void command(int, char **);

 private:
  int *dlist;
  int compress_flag,mol_flag;
  std::map<tagint,int> *hash;

  void delete_group(int, char **);
  void delete_region(int, char **);
  void delete_overlap(int, char **);
  void delete_porosity(int, char **);
  void delete_molecule();
  void recount_topology();
  void options(int, char **);

  inline int sbmask(int j) {
    return j >> SBBITS & 3;
  }

  // static variable for ring communication callback to access class data
  // callback functions for ring communication

  static DeleteAtoms *cptr;
  static void molring(int, char *);
};

}

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

E: Could not find delete_atoms group ID

Group ID used in the delete_atoms command does not exist.

E: Could not find delete_atoms region ID

Region ID used in the delete_atoms command does not exist.

E: Delete_atoms requires a pair style be defined

This is because atom deletion within a cutoff uses a pairwise
neighbor list.

E: Delete_atoms cutoff > neighbor cutoff

Cannot delete atoms further away than a processor knows about.

*/

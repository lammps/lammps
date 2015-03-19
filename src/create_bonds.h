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

CommandStyle(create_bonds,CreateBonds)

#else

#ifndef LMP_CREATE_BONDS_H
#define LMP_CREATE_BONDS_H

#include "pointers.h"

namespace LAMMPS_NS {

class CreateBonds : protected Pointers {
 public:
  CreateBonds(class LAMMPS *);
  void command(int, char **);

 private:
  inline int sbmask(int j) {
    return j >> SBBITS & 3;
  }
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Create_bonds command before simulation box is defined

UNDOCUMENTED

E: Cannot use create_bonds unless atoms have IDs

UNDOCUMENTED

E: Cannot use create_bonds with non-molecular system

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot find create_bonds group ID

UNDOCUMENTED

E: Invalid bond type in create_bonds command

UNDOCUMENTED

E: Create_bonds requires a pair style be defined

UNDOCUMENTED

E: Create_bonds max distance > neighbor cutoff

UNDOCUMENTED

W: Create_bonds max distance > minimum neighbor cutoff

UNDOCUMENTED

E: Create_bonds command requires special_bonds 1-2 weights be 0.0

UNDOCUMENTED

E: Create_bonds command requires no kspace_style be defined

UNDOCUMENTED

E: New bond exceeded bonds per atom in create_bonds

UNDOCUMENTED

U: Delete_atoms command before simulation box is defined

The delete_atoms command cannot be used before a read_data,
read_restart, or create_box command.

U: Cannot use delete_atoms unless atoms have IDs

Your atoms do not have IDs, so the delete_atoms command cannot be
used.

U: Could not find delete_atoms group ID

Group ID used in the delete_atoms command does not exist.

U: Could not find delete_atoms region ID

Region ID used in the delete_atoms command does not exist.

U: Delete_atoms requires a pair style be defined

This is because atom deletion within a cutoff uses a pairwise
neighbor list.

U: Delete_atoms cutoff > neighbor cutoff

Cannot delete atoms further away than a processor knows about.

*/

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
  int igroup,group1bit,group2bit;
  int btype,atype,dtype;
  tagint batom1,batom2,aatom1,aatom2,aatom3,datom1,datom2,datom3,datom4;
  double rmin,rmax;

  void many();
  void single_bond();
  void single_angle();
  void single_dihedral();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Create_bonds command before simulation box is defined

Self-explanatory.

E: Cannot use create_bonds unless atoms have IDs

This command requires a mapping from global atom IDs to local atoms,
but the atoms that have been defined have no IDs.

E: Cannot use create_bonds with non-molecular system

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot find create_bonds group ID

Self-explanatory.

E: Invalid bond type in create_bonds command

Self-explanatory.

E: Cannot use special no with create_bonds many

UNDOCUMENTED

E: Invalid angle type in create_bonds command

UNDOCUMENTED

E: Invalid dihedral type in create_bonds command

UNDOCUMENTED

E: Create_bonds requires a pair style be defined

Self-explanatory.

E: Create_bonds max distance > neighbor cutoff

Can only create bonds for atom pairs that will be in neighbor list.

W: Create_bonds max distance > minimum neighbor cutoff

This means atom pairs for some atom types may not be in the neighbor
list and thus no bond can be created between them.

E: Create_bonds command requires special_bonds 1-2 weights be 0.0

This is so that atom pairs that are already bonded to not appear in
the neighbor list.

E: Create_bonds command requires no kspace_style be defined

This is so that atom pairs that are already bonded to not appear
in the neighbor list.

E: New bond exceeded bonds per atom in create_bonds

See the read_data command for info on setting the "extra bond per
atom" header value to allow for additional bonds to be formed.

E: Create_bonds single/bond atoms do not exist

UNDOCUMENTED

E: Create_bonds single/angle atoms do not exist

UNDOCUMENTED

E: New angle exceeded angles per atom in create_bonds

UNDOCUMENTED

E: Create_bonds single/dihedral atoms do not exist

UNDOCUMENTED

E: New dihedral exceeded dihedrals per atom in create_bonds

UNDOCUMENTED

*/

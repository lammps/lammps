/* ----------------------------------------------------------------------
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

CommandStyle(create_atoms,CreateAtoms)

#else

#ifndef LMP_CREATE_ATOMS_H
#define LMP_CREATE_ATOMS_H

#include "pointers.h"

namespace LAMMPS_NS {

class CreateAtoms : protected Pointers {
 public:
  CreateAtoms(class LAMMPS *);
  void command(int, char **);

 private:
  int ntype,style,mode,nregion,nbasis,nrandom,seed;
  int *basistype;
  double xone[3];
  int remapflag;

  class Molecule *onemol;
  class RanMars *ranmol;

  int triclinic;
  double sublo[3],subhi[3];   // epsilon-extended proc sub-box for adding atoms

  void add_single();
  void add_random();
  void add_lattice();
  void add_molecule(double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Create_atoms command before simulation box is defined

The create_atoms command cannot be used before a read_data,
read_restart, or create_box command.

E: Cannot create_atoms after reading restart file with per-atom info

The per-atom info was stored to be used when by a fix that you
may re-define.  If you add atoms before re-defining the fix, then
there will not be a correct amount of per-atom info.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid atom type in create_atoms command

The create_box command specified the range of valid atom types.
An invalid type is being requested.

E: Create_atoms region ID does not exist

A region ID used in the create_atoms command does not exist.

E: Invalid basis setting in create_atoms command

UNDOCUMENTED

E: Cannot create atoms with undefined lattice

Must use the lattice command before using the create_atoms
command.

E: Too many total atoms

See the setting for bigint in the src/lmptype.h file.

E: No overlap of box and region for create_atoms

Self-explanatory.

*/

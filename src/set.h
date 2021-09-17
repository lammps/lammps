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
CommandStyle(set,Set);
// clang-format on
#else

#ifndef LMP_SET_H
#define LMP_SET_H

#include "command.h"

namespace LAMMPS_NS {

class Set : public Command {
 public:
  Set(class LAMMPS *lmp) : Command(lmp){};
  void command(int, char **);

 private:
  char *id;
  int *select;
  int style, ivalue, newtype, count, index_custom, icol_custom;
  int ximage, yimage, zimage, ximageflag, yimageflag, zimageflag;
  int cc_index;
  bigint nsubset;
  double dvalue, xvalue, yvalue, zvalue, wvalue, fraction;

  int varflag, varflag1, varflag2, varflag3, varflag4;
  int ivar1, ivar2, ivar3, ivar4;
  double *vec1, *vec2, *vec3, *vec4;

  int discflag;

  void selection(int);
  void set(int);
  void setrandom(int);
  void topology(int);
  void varparse(const char *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Set command before simulation box is defined

The set command cannot be used before a read_data, read_restart,
or create_box command.

E: Set command with no atoms existing

No atoms are yet defined so the set command cannot be used.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid value in set command

The value specified for the setting is invalid, likely because it is
too small or too large.

E: Invalid random number seed in set command

Random number seed must be > 0.

E: Cannot set this attribute for this atom style

The attribute being set does not exist for the defined atom style.

E: Invalid dipole length in set command

Self-explanatory.

E: Invalid density in set command

Density must be > 0.0.

E: Density/disc option requires 2d simulation

UNDOCUMENTED

E: Invalid volume in set command

Volume must be > 0.0.

E: Cannot set non-zero image flag for non-periodic dimension

Self-explanatory.

E: Cannot set meso/e for this atom style

Self-explanatory.

E: Cannot set meso/cv for this atom style

Self-explanatory.

E: Cannot set meso/rho for this atom style

Self-explanatory.

E: Cannot set edpd/temp for this atom style

UNDOCUMENTED

E: Cannot set edpd/cv for this atom style

UNDOCUMENTED

E: Cannot set cc for this atom style

UNDOCUMENTED

E: Cannot set smd/mass/density for this atom style

Self-explanatory.

E: Cannot set smd/contact/radius for this atom style

Self-explanatory.

E: Cannot set dpd/theta for this atom style

Self-explanatory.

E: Set command integer vector does not exist

Self-explanatory.

E: Set command floating point vector does not exist

Self-explanatory.

E: Cannot use set atom with no atom IDs defined

Atom IDs are not defined, so they cannot be used to identify an atom.

E: Cannot use set mol with no molecule IDs defined

Self-explanatory.

E: Could not find set group ID

Group ID specified in set command does not exist.

E: Set region ID does not exist

Region ID specified in set command does not exist.

W: Changing a property of atoms in rigid bodies that has no effect unless rigid bodies are rebuild

UNDOCUMENTED

E: Invalid mass in set command

Self-explanatory.

E: Invalid diameter in set command

Self-explanatory.

E: Invalid shape in set command

Self-explanatory.

E: Invalid length in set command

Self-explanatory.

E: Cannot set quaternion for atom that has none

Self-explanatory.

E: Cannot set quaternion with xy components for 2d system

Self-explanatory.

E: Cannot set theta for atom that is not a line

Self-explanatory.

E: Cannot set bond topology types for atom style template

The bond, angle, etc types cannot be changed for this atom style since
they are static settings in the molecule template files.

E: Bond atom missing in set command

The set command cannot find one or more atoms in a particular bond on
a particular processor.  The pairwise cutoff is too short or the atoms
are too far apart to make a valid bond.

E: Angle atom missing in set command

The set command cannot find one or more atoms in a particular angle on
a particular processor.  The pairwise cutoff is too short or the atoms
are too far apart to make a valid angle.

E: Dihedral atom missing in set command

The set command cannot find one or more atoms in a particular dihedral
on a particular processor.  The pairwise cutoff is too short or the
atoms are too far apart to make a valid dihedral.

E: Improper atom missing in set command

The set command cannot find one or more atoms in a particular improper
on a particular processor.  The pairwise cutoff is too short or the
atoms are too far apart to make a valid improper.

E: Variable name for set command does not exist

Self-explanatory.

E: Variable for set command is invalid style

Only atom-style variables can be used.

*/

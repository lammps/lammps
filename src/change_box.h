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

CommandStyle(change_box,ChangeBox)

#else

#ifndef LMP_CHANGE_BOX_H
#define LMP_CHANGE_BOX_H

#include "pointers.h"

namespace LAMMPS_NS {

class ChangeBox : protected Pointers {
 public:
  ChangeBox(class LAMMPS *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Change_box command before simulation box is defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Change_box operation is invalid

Cannot change orthogonal box to orthogonal or a triclinic box to
triclinic.

E: Cannot change box to orthogonal when tilt is non-zero

Self-explanatory

E: Cannot change box with dumps defined

Self-explanatory.

E: Cannot change box with certain fixes defined

The change_box command cannot be used when fix ave/spatial or
fix/deform are defined .

*/

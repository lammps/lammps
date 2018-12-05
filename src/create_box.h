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

CommandStyle(create_box,CreateBox)

#else

#ifndef LMP_CREATE_BOX_H
#define LMP_CREATE_BOX_H

#include "pointers.h"

namespace LAMMPS_NS {

class CreateBox : protected Pointers {
 public:
  CreateBox(class LAMMPS *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot create_box after simulation box is defined

A simulation box can only be defined once.

E: Cannot run 2d simulation with non-periodic Z dimension

Use the boundary command to make the z dimension periodic in order to
run a 2d simulation.

E: Create_box region ID does not exist

Self-explanatory.

E: Create_box region does not support a bounding box

Not all regions represent bounded volumes.  You cannot use
such a region with the create_box command.

E: No bonds allowed with this atom style

Self-explanatory.

E: No angles allowed with this atom style

Self-explanatory.

E: No dihedrals allowed with this atom style

Self-explanatory.

E: No impropers allowed with this atom style

Self-explanatory.

*/

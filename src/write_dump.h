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

CommandStyle(write_dump,WriteDump)

#else

#ifndef LMP_WRITE_DUMP_H
#define LMP_WRITE_DUMP_H

#include "pointers.h"

namespace LAMMPS_NS {

class WriteDump : protected Pointers {
 public:
  WriteDump(class LAMMPS *lmp) : Pointers(lmp) {};
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

E: Unrecognized dump style

The choice of dump style is unknown.

W: Calling write_dump before a full system init.

The write_dump command is used before the system has been fully
initialized as part of a 'run' or 'minimize' command. Not all dump
styles and features are fully supported at this point and thus the
command may fail or produce incomplete or incorrect output. Insert
a "run 0" command, if a full system init is required.

*/

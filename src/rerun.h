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

CommandStyle(rerun,Rerun)

#else

#ifndef LMP_RERUN_H
#define LMP_RERUN_H

#include "pointers.h"

namespace LAMMPS_NS {

class Rerun : protected Pointers {
 public:
  Rerun(class LAMMPS *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Rerun command before simulation box is defined

The rerun command cannot be used before a read_data, read_restart, or
create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Rerun dump file does not contain requested snapshot

Self-explanatory.

*/

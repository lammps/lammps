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

CommandStyle(run,Run);

#else

#ifndef LMP_RUN_H
#define LMP_RUN_H

#include "pointers.h"

namespace LAMMPS_NS {

class Run : protected Pointers {
 public:
  Run(class LAMMPS *);
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

E: Run command before simulation box is defined

The run command cannot be used before a read_data, read_restart, or
create_box command.

E: Invalid run command N value

The number of timesteps must fit in a 32-bit integer.  If you want to
run for more steps than this, perform multiple shorter runs.

E: Invalid run command upto value

Self-explanatory.

E: Invalid run command start/stop value

Self-explanatory.

E: Run command start value is after start of run

Self-explanatory.

E: Run command stop value is before end of run

Self-explanatory.

E: Run flag 'pre no' not compatible with r-RESPA

UNDOCUMENTED

E: Too many timesteps

The cumulative timesteps must fit in a 64-bit integer.

*/

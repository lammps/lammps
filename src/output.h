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

#ifndef LMP_OUTPUT_H
#define LMP_OUTPUT_H

#include "pointers.h"

namespace LAMMPS_NS {

class Output : protected Pointers {
 public:
  bigint next;                 // next timestep for any kind of output

  bigint next_thermo;          // next timestep for thermo output
  int thermo_every;            // thermo output every this many steps
  bigint last_thermo;          // last timestep thermo was output
  char *var_thermo;            // variable name for thermo frequency
  int ivar_thermo;             // variable index for thermo frequency
  class Thermo *thermo;        // Thermodynamic computations

  int ndump;                   // # of Dumps defined
  int max_dump;                // max size of Dump list
  bigint next_dump_any;        // next timestep for any Dump
  int *every_dump;             // output of each Dump every this many steps
  bigint *next_dump;           // next timestep to do each Dump
  bigint *last_dump;           // last timestep each snapshot was output
  char **var_dump;             // variable name for dump frequency
  int *ivar_dump;              // variable index for dump frequency
  class Dump **dump;           // list of defined Dumps

  bigint next_restart;         // next timestep to write a restart file
  int restart_every;           // write a restart file every this many steps
  bigint last_restart;         // last timestep a restart file was output
  int restart_toggle;          // 0 if use restart1 as prefix
                               // 1 if use restart1 as file, 2 for restart2
  char *restart1,*restart2;    // names for restart files
  class WriteRestart *restart; // Restart output

  Output(class LAMMPS *);
  ~Output();
  void init();
  void setup(int);                   // initial output before run/min
  void write(bigint);                // output for current timestep
  void write_dump(bigint);           // force output of dump snapshots
  void write_restart(bigint);        // force output of a restart file

  void add_dump(int, char **);       // add a Dump to Dump list
  void modify_dump(int, char **);    // modify a Dump
  void delete_dump(char *);          // delete a Dump from Dump list

  void create_thermo(int, char **);  // create a thermo style
  void create_restart(int, char **); // create Restart and restart files

  void memory_usage();               // print out memory usage
};

}

#endif

/* ERROR/WARNING messages:

E: Variable name for thermo every does not exist

Self-explanatory.

E: Variable for thermo every is invalid style

Only equal-style variables can be used.

E: Variable name for dump every does not exist

Self-explanatory.

E: Variable for dump every is invalid style

Only equal-style variables can be used.

E: Dump every variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

E: Thermo every variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Reuse of dump ID

A dump ID cannot be used twice.

E: Could not find dump group ID

A group ID used in the dump command does not exist.

E: Invalid dump frequency

Dump frequency must be 1 or greater.

E: Invalid dump style

The choice of dump style is unknown.

E: Cound not find dump_modify ID

Self-explanatory.

E: Could not find undump ID

A dump ID used in the undump command does not exist.

E: Thermo_style command before simulation box is defined

The thermo_style command cannot be used before a read_data,
read_restart, or create_box command.

W: New thermo_style command, previous thermo_modify settings will be lost

If a thermo_style command is used after a thermo_modify command, the
settings changed by the thermo_modify command will be reset to their
default values.  This is because the thermo_modify commmand acts on
the currently defined thermo style, and a thermo_style command creates
a new style.

*/

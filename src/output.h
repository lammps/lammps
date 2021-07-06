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

#ifndef LMP_OUTPUT_H
#define LMP_OUTPUT_H

#include "pointers.h"

#include <map>

namespace LAMMPS_NS {

class Dump;

class Output : protected Pointers {
 public:
  bigint next;    // next timestep for any kind of output

  bigint next_thermo;      // next timestep for thermo output
  int thermo_every;        // output freq for thermo, 0 if first/last only
  bigint last_thermo;      // last timestep thermo was output
  char *var_thermo;        // variable name for thermo freq, null pointer if every
  int ivar_thermo;         // variable index for thermo frequency
  class Thermo *thermo;    // Thermodynamic computations

  int ndump;               // # of Dumps defined
  int max_dump;            // max size of Dump list
  bigint next_dump_any;    // next timestep for any Dump
  int *every_dump;         // write freq for each Dump, 0 if var
  bigint *next_dump;       // next timestep to do each Dump
  bigint *last_dump;       // last timestep each snapshot was output
  char **var_dump;         // variable name for dump frequency
  int *ivar_dump;          // variable index for dump frequency
  Dump **dump;             // list of defined Dumps

  int restart_flag;               // 1 if any restart files are written
  int restart_flag_single;        // 1 if single restart files are written
  int restart_flag_double;        // 1 if double restart files are written
  bigint next_restart;            // next timestep to write any restart file
  bigint next_restart_single;     // next timestep to write a single restart file
  bigint next_restart_double;     // next timestep to write a double restart file
  int restart_every_single;       // single restart file write freq, 0 if var
  int restart_every_double;       // double restart file write freq, 0 if var
  bigint last_restart;            // last timestep any restart file was output
  int restart_toggle;             // 0 if use restart2a as prefix, 1 if restart2b
  char *var_restart_single;       // variable name for single restart freq
  char *var_restart_double;       // variable name for double restart freq
  int ivar_restart_single;        // index of var_restart_single
  int ivar_restart_double;        // index of var_restart_double
  char *restart1;                 // name single restart file
  char *restart2a, *restart2b;    // names of double restart files
  class WriteRestart *restart;    // class for writing restart files

  typedef Dump *(*DumpCreator)(LAMMPS *, int, char **);
  typedef std::map<std::string, DumpCreator> DumpCreatorMap;
  DumpCreatorMap *dump_map;

  Output(class LAMMPS *);
  ~Output();
  void init();
  void setup(int memflag = 1);    // initial output before run/min
  void write(bigint);             // output for current timestep
  void write_dump(bigint);        // force output of dump snapshots
  void write_restart(bigint);     // force output of a restart file
  void reset_timestep(bigint);    // reset next timestep for all output

  void add_dump(int, char **);       // add a Dump to Dump list
  void modify_dump(int, char **);    // modify a Dump
  void delete_dump(char *);          // delete a Dump from Dump list
  int find_dump(const char *);       // find a Dump ID

  void set_thermo(int, char **);        // set thermo output freqquency
  void create_thermo(int, char **);     // create a thermo style
  void create_restart(int, char **);    // create Restart and restart files

  void memory_usage();    // print out memory usage

 private:
  template <typename T> static Dump *dump_creator(LAMMPS *, int, char **);
};

}    // namespace LAMMPS_NS

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

E: Variable name for restart does not exist

Self-explanatory.

E: Variable for restart is invalid style

Only equal-style variables can be used.

E: Dump every variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

E: Restart variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

E: Thermo every variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

E: Thermo_modify every variable returned a bad timestep

The returned timestep is less than or equal to the current timestep.

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

E: Unrecognized dump style

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
default values.  This is because the thermo_modify command acts on
the currently defined thermo style, and a thermo_style command creates
a new style.

E: Both restart files must use % or neither

Self-explanatory.

E: Both restart files must use MPI-IO or neither

Self-explanatory.

*/

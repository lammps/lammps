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

#ifndef OUTPUT_H
#define OUTPUT_H

#include "pointers.h"

namespace LAMMPS_NS {

class Output : protected Pointers {
 public:
  int next;                    // next timestep for any kind of output

  int next_thermo;             // next timestep for thermo output
  int thermo_every;            // thermo output every this many steps
  int last_thermo;             // last timestep thermo was output
  class Thermo *thermo;        // Thermodynamic computations

  int ndump;                   // # of Dumps defined
  int max_dump;                // max size of Dump list
  int next_dump_any;           // next timestep for any Dump
  int *next_dump;              // next timestep to do each Dump
  int *dump_every;             // output of each Dump every this many steps
  int *last_dump;              // last timestep each a snapshot was output
  class Dump **dump;           // list of defined Dumps

  int next_restart;            // next timestep to write a restart file
  int restart_every;           // write a restart file every this many steps
  int last_restart;            // last timestep a restart file was output
  int restart_toggle;          // 0 if use restart1 as prefix
                               // 1 if use restart1 as file, 2 for restart2
  char *restart1,*restart2;    // names for restart files
  class WriteRestart *restart; // Restart output

  Output(class LAMMPS *);
  ~Output();
  void init();
  void setup(int);                   // initial output before run/min
  void write(int);                   // output for current timestep
  void write_dump(int);              // force output of dump snapshots
  void write_restart(int);           // force output of a restart file

  void add_dump(int, char **);       // add a Dump to Dump list
  void modify_dump(int, char **);    // modify a Dump
  void delete_dump(char *);          // delete a Dump from Dump list

  void create_thermo(int, char **);  // create a thermo style
  void create_restart(int, char **); // create Restart and restart files

  void memory_usage();               // print out memory usage
};

}

#endif

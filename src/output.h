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

  int ndump;                    // # of Dumps defined
  int max_dump;                 // max size of Dump list
  bigint next_dump_any;         // next timestep for any dump
  bigint next_time_dump_any;    // next timestep for any time dump with computes
  int any_time_dumps;           // 1 if any time dump defined
  int *mode_dump;               // 0/1 if write every N timesteps or Delta in sim time
  int *every_dump;              // dump every N timesteps, 0 if variable
  double *every_time_dump;      // dump every Delta of sim time, 0.0 if variable
  bigint *next_dump;            // next timestep to perform dump
  double *next_time_dump;       // next simulation time to perform dump (mode = 1)
  bigint *last_dump;            // last timestep each snapshot was output
  char **var_dump;              // variable name for next dump (steps or sim time)
  int *ivar_dump;               // variable index of var_dump name
  Dump **dump;                  // list of defined Dumps

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
  ~Output() override;
  void init();
  void setup(int memflag = 1);    // initial output before run/min
  void write(bigint);             // output for current timestep
  void write_dump(bigint);        // force output of dump snapshots
  void write_restart(bigint);     // force output of a restart file
  void reset_timestep(bigint);    // reset output which depends on timestep
  void reset_dt();                // reset output which depends on timestep size

  Dump *add_dump(int, char **);                       // add a Dump to Dump list
  void modify_dump(int, char **);                     // modify a Dump
  void delete_dump(const std::string &);              // delete a Dump from Dump list
  Dump *get_dump_by_id(const std::string &) const;    // find a Dump by ID
  Dump *get_dump_by_index(int idx) const              // find a Dump by index in Dump list
  {
    return ((idx >= 0) && (idx < ndump)) ? dump[idx] : nullptr;
  }

  const std::vector<Dump *> &get_dump_list();    // get vector with all dumps
  int check_time_dumps(bigint);                  // check if any time dump is output now

  void set_thermo(int, char **);        // set thermo output freqquency
  void create_thermo(int, char **);     // create a thermo style
  void create_restart(int, char **);    // create Restart and restart files

  void memory_usage();    // print out memory usage

 private:
  std::vector<Dump *> dump_list;
  void calculate_next_dump(int, int, bigint);
};

}    // namespace LAMMPS_NS

#endif

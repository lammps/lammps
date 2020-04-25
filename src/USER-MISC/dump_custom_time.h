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

#ifdef DUMP_CLASS

DumpStyle(custom/time,DumpCustom_Time)

#else

#ifndef LMP_DUMP_CUSTOM_TIME_H
#define LMP_DUMP_CUSTOM_TIME_H

#include "dump_custom.h"

namespace LAMMPS_NS {

class DumpCustom_Time : public DumpCustom {
 public:
  DumpCustom_Time(class LAMMPS *, int, char **);
  virtual ~DumpCustom_Time();

 protected:
  double time_every; //time frequency of a dump
  double next_time;  //next time for a dump
  double tol;
  int write_time; // -1 - write in the next timestep, 0 - write, >0 check if time for write

  // private methods

//  virtual void init_style();
  virtual void write_header(bigint);
  virtual void write_data(int, double *);
  int count();
//  void pack(tagint *);
//  virtual int convert_string(int, double *);
//  bigint memory_usage();
//
//  int parse_fields(int, char **);
//  int add_compute(char *);
//  int add_fix(char *);
//  int add_variable(char *);
//  int add_custom(char *, int);
//  virtual int modify_param(int, char **);
//
//  typedef void (DumpCustom_Time::*FnPtrHeader)(bigint);
//  FnPtrHeader header_choice;           // ptr to write header functions
//  void header_binary(bigint);
//  void header_binary_triclinic(bigint);
//  void header_item(bigint);
//  void header_item_triclinic(bigint);
//
//  typedef int (DumpCustom_Time::*FnPtrConvert)(int, double *);
//  FnPtrConvert convert_choice;          // ptr to convert data functions
//  int convert_image(int, double *);
//  int convert_noimage(int, double *);
//
//  typedef void (DumpCustom_Time::*FnPtrWrite)(int, double *);
//  FnPtrWrite write_choice;             // ptr to write data functions
//  void write_binary(int, double *);
//  void write_string(int, double *);
//  void write_lines(int, double *);
//
//  // customize by adding a method prototype
//
//  typedef void (DumpCustom_Time::*FnPtrPack)(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: No dump custom arguments specified

The dump custom command requires that atom quantities be specified to
output to dump file.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid attribute in dump custom command

Self-explanatory.

E: Dump_modify format line is too short

UNDOCUMENTED

E: Could not find dump custom compute ID

Self-explanatory.

E: Could not find dump custom fix ID

Self-explanatory.

E: Dump custom and fix not computed at compatible times

The fix must produce per-atom quantities on timesteps that dump custom
needs them.

E: Could not find dump custom variable name

Self-explanatory.

E: Could not find custom per-atom property ID

Self-explanatory.

E: Region ID for dump custom does not exist

Self-explanatory.

E: Compute used in dump between runs is not current

The compute was not invoked on the current timestep, therefore it
cannot be used in a dump between runs.

E: Threshold for an atom property that isn't allocated

A dump threshold has been requested on a quantity that is
not defined by the atom style used in this simulation.

E: Dumping an atom property that isn't allocated

The chosen atom style does not define the per-atom quantity being
dumped.

E: Dump custom compute does not compute per-atom info

Self-explanatory.

E: Dump custom compute does not calculate per-atom vector

Self-explanatory.

E: Dump custom compute does not calculate per-atom array

Self-explanatory.

E: Dump custom compute vector is accessed out-of-range

Self-explanatory.

E: Dump custom fix does not compute per-atom info

Self-explanatory.

E: Dump custom fix does not compute per-atom vector

Self-explanatory.

E: Dump custom fix does not compute per-atom array

Self-explanatory.

E: Dump custom fix vector is accessed out-of-range

Self-explanatory.

E: Dump custom variable is not atom-style variable

Only atom-style variables generate per-atom quantities, needed for
dump output.

E: Custom per-atom property ID is not floating point

Self-explanatory.

E: Custom per-atom property ID is not integer

Self-explanatory.

E: Dump_modify region ID does not exist

Self-explanatory.

E: Dump_modify int format does not contain d character

UNDOCUMENTED

E: Dump_modify element names do not match atom types

UNDOCUMENTED

E: Dump modify can only have one refresh

UNDOCUMENTED

E: Invalid attribute in dump modify command

Self-explanatory.

E: Could not find dump modify compute ID

Self-explanatory.

E: Dump modify compute ID does not compute per-atom info

Self-explanatory.

E: Dump modify compute ID does not compute per-atom vector

Self-explanatory.

E: Dump modify compute ID does not compute per-atom array

Self-explanatory.

E: Dump modify compute ID vector is not large enough

Self-explanatory.

E: Could not find dump modify fix ID

Self-explanatory.

E: Dump modify fix ID does not compute per-atom info

Self-explanatory.

E: Dump modify fix ID does not compute per-atom vector

Self-explanatory.

E: Dump modify fix ID does not compute per-atom array

Self-explanatory.

E: Dump modify fix ID vector is not large enough

Self-explanatory.

E: Could not find dump modify variable name

Self-explanatory.

E: Dump modify variable is not atom-style variable

Self-explanatory.

E: Could not find dump modify custom atom floating point property ID

Self-explanatory.

E: Could not find dump modify custom atom integer property ID

Self-explanatory.

E: Invalid dump_modify thresh attribute

UNDOCUMENTED

E: Invalid dump_modify thresh operator

UNDOCUMENTED

U: Dump_modify format string is too short

There are more fields to be dumped in a line of output than your
format string specifies.

U: Dump modify element names do not match atom types

Number of element names must equal number of atom types.

U: Invalid dump_modify threshold operator

Operator keyword used for threshold specification in not recognized.

*/

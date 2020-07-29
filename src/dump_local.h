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

DumpStyle(local,DumpLocal)

#else

#ifndef LMP_DUMP_LOCAL_H
#define LMP_DUMP_LOCAL_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpLocal : public Dump {
 public:
  DumpLocal(LAMMPS *, int, char **);
  virtual ~DumpLocal();

 protected:
  int nevery;                // dump frequency to check Fix against
  char *label;               // string for dump file header

  int nmine;                 // # of lines I am dumping
  int *vtype;                // type of each vector (INT, DOUBLE)
  char **vformat;            // format string for each vector element

  char *columns;             // column labels

  int nfield;                // # of keywords listed by user

  int *field2index;          // which compute,fix,variable calcs this field
  int *argindex;             // index into compute,fix scalar_atom,vector_atom
                             // 0 for scalar_atom, 1-N for vector_atom values

  int ncompute;              // # of Compute objects used by dump
  char **id_compute;         // their IDs
  class Compute **compute;   // list of ptrs to the Compute objects

  int nfix;                  // # of Fix objects used by dump
  char **id_fix;             // their IDs
  class Fix **fix;           // list of ptrs to the Fix objects

  void init_style();
  int modify_param(int, char **);
  virtual void write_header(bigint);
  int count();
  void pack(tagint *);
  int convert_string(int, double *);
  virtual void write_data(int, double *);

  void parse_fields(int, char **);
  int add_compute(char *);
  int add_fix(char *);

  typedef void (DumpLocal::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;             // ptr to write data functions
  void write_string(int, double *);
  void write_lines(int, double *);

  // customize by adding a method prototype

  typedef void (DumpLocal::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_index(int);
  void pack_compute(int);
  void pack_fix(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: No dump local arguments specified

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Binary files are not supported with dump local

UNDOCUMENTED

E: Dump local cannot sort by atom ID

This is because dump local does not really dump per-atom info.

E: Dump_modify format line is too short

UNDOCUMENTED

E: Could not find dump local compute ID

Self-explanatory.

E: Could not find dump local fix ID

Self-explanatory.

E: Dump local and fix not computed at compatible times

The fix must produce per-atom quantities on timesteps that dump local
needs them.

E: Compute used in dump between runs is not current

The compute was not invoked on the current timestep, therefore it
cannot be used in a dump between runs.

E: Dump local count is not consistent across input fields

Every column of output must be the same length.

E: Invalid attribute in dump local command

Self-explanatory.

E: Dump local compute does not compute local info

Self-explanatory.

E: Dump local compute does not calculate local vector

Self-explanatory.

E: Dump local compute does not calculate local array

Self-explanatory.

E: Dump local compute vector is accessed out-of-range

Self-explanatory.

E: Dump local fix does not compute local info

Self-explanatory.

E: Dump local fix does not compute local vector

Self-explanatory.

E: Dump local fix does not compute local array

Self-explanatory.

E: Dump local fix vector is accessed out-of-range

Self-explanatory.

E: Dump local attributes contain no compute or fix

Self-explanatory.

*/

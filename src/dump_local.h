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
  ~DumpLocal();
  void init();

 private:
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

  int modify_param(int, char **);
  void write_header(int);
  int count();
  int pack();
  void write_data(int, double *);

  void parse_fields(int, char **);
  int add_compute(char *);
  int add_fix(char *);

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

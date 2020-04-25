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

DumpStyle(xyz/time,DumpXYZ_Time)

#else

#ifndef LMP_DUMP_XYZ_TIME_H
#define LMP_DUMP_XYZ_TIME_H

#include "dump_xyz.h"

#include <cstring>
#include <cstdio>

char** RemoveTimeIncrement(int narg, char**  arg);

namespace LAMMPS_NS {

class DumpXYZ_Time : public DumpXYZ {
 public:
  DumpXYZ_Time(class LAMMPS *, int, char**);
  virtual ~DumpXYZ_Time();

 protected:
  double time_every; //time frequency of a dump
  double next_time;  //next time for a dump
  double tol;    
  bool write_time;

  void write_header(bigint);
  void write_data(int, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid dump xyz filename

Filenames used with the dump xyz style cannot be binary or cause files
to be written by each processor.

E: Dump modify element names do not match atom types

Number of element names must equal number of atom types.

*/

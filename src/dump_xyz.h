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

DumpStyle(xyz,DumpXYZ)

#else

#ifndef LMP_DUMP_XYZ_H
#define LMP_DUMP_XYZ_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpXYZ : public Dump {
 public:
  DumpXYZ(class LAMMPS *, int, char**);
  ~DumpXYZ() {}

 private:
  void init_style();
  void write_header(bigint);
  int count();
  void pack(int *);
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

*/

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

#ifdef FIX_CLASS

FixStyle(vector,FixVector)

#else

#ifndef LMP_FIX_VECTOR_H
#define LMP_FIX_VECTOR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixVector : public Fix {
 public:
  FixVector(class LAMMPS *, int, char **);
  ~FixVector();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_vector(int);
  double compute_array(int,int);

 private:
  int nvalues;
  int *which,*argindex,*value2index;
  char **ids;

  bigint nextstep,initialstep;

  int ncount;        // # of values currently in growing vector or array
  int ncountmax;     // max # of values vector/array can hold
  double *vector;
  double **array;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Compute ID for fix vector does not exist

UNDOCUMENTED

E: Fix vector compute does not calculate a scalar

UNDOCUMENTED

E: Fix vector compute does not calculate a vector

UNDOCUMENTED

E: Fix vector compute vector is accessed out-of-range

UNDOCUMENTED

E: Fix ID for fix vector does not exist

UNDOCUMENTED

E: Fix vector fix does not calculate a scalar

UNDOCUMENTED

E: Fix vector fix does not calculate a vector

UNDOCUMENTED

E: Fix vector fix vector is accessed out-of-range

UNDOCUMENTED

E: Fix for fix vector not computed at compatible time

UNDOCUMENTED

E: Variable name for fix vector does not exist

UNDOCUMENTED

E: Fix vector variable is not equal-style variable

UNDOCUMENTED

E: Fix vector cannot set output array intensive/extensive from these inputs

UNDOCUMENTED

E: Overflow of allocated fix vector storage

UNDOCUMENTED

*/

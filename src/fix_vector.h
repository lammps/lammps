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

#ifdef FIX_CLASS
// clang-format off
FixStyle(vector,FixVector);
// clang-format on
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
  double compute_array(int, int);

 private:
  int nvalues;
  int *which, *argindex, *value2index;
  char **ids;

  bigint nextstep, initialstep;

  int ncount;       // # of values currently in growing vector or array
  int ncountmax;    // max # of values vector/array can hold
  double *vector;
  double **array;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for fix vector does not exist

Self-explanatory.

E: Fix vector compute does not calculate a scalar

Self-explanatory.

E: Fix vector compute does not calculate a vector

Self-explanatory.

E: Fix vector compute vector is accessed out-of-range

Self-explanatory.

E: Fix ID for fix vector does not exist

Self-explanatory.

E: Fix vector fix does not calculate a scalar

Self-explanatory.

E: Fix vector fix does not calculate a vector

Self-explanatory.

E: Fix vector fix vector is accessed out-of-range

Self-explanatory.

E: Fix for fix vector not computed at compatible time

Fixes generate their values on specific timesteps.  Fix vector is
requesting a value on a non-allowed timestep.

E: Variable name for fix vector does not exist

Self-explanatory.

E: Fix vector variable is not equal-style variable

Self-explanatory.

E: Fix vector variable is not vector-style variable

UNDOCUMENTED

E: Fix vector cannot set output array intensive/extensive from these inputs

The inputs to the command have conflicting intensive/extensive attributes.
You need to use more than one fix vector command.

E: Overflow of allocated fix vector storage

This should not normally happen if the fix correctly calculated
how long the vector will grow to.  Contact the developers.

*/

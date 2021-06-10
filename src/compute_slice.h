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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(slice,ComputeSlice);
// clang-format on
#else

#ifndef LMP_COMPUTE_SLICE_H
#define LMP_COMPUTE_SLICE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSlice : public Compute {
 public:
  ComputeSlice(class LAMMPS *, int, char **);
  virtual ~ComputeSlice();
  void init();
  void compute_vector();
  void compute_array();

 private:
  int me;
  int nstart, nstop, nskip, nvalues;
  int *which, *argindex, *value2index;
  char **ids;

  void extract_one(int, double *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for compute slice does not exist

Self-explanatory.

E: Compute slice compute does not calculate a global array

Self-explanatory.

E: Compute slice compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Compute slice compute does not calculate a global vector

Self-explanatory.

E: Compute slice compute array is accessed out-of-range

An index for the array is out of bounds.

E: Compute slice compute does not calculate global vector or array

Self-explanatory.

E: Fix ID for compute slice does not exist

Self-explanatory.

E: Compute slice fix does not calculate a global array

Self-explanatory.

E: Compute slice fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Compute slice fix does not calculate a global vector

Self-explanatory.

E: Compute slice fix array is accessed out-of-range

An index for the array is out of bounds.

E: Compute slice fix does not calculate global vector or array

Self-explanatory.

E: Variable name for compute slice does not exist

UNDOCUMENTED

E: Compute slice variable is not vector-style variable

UNDOCUMENTED

E: Compute slice vector variable cannot be indexed

UNDOCUMENTED

E: Fix used in compute slice not computed at compatible time

Fixes generate their values on specific timesteps.  Compute slice is
requesting a value on a non-allowed timestep.

E: Compute slice variable is not long enough

UNDOCUMENTED

*/

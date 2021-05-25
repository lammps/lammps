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
FixStyle(ave/correlate,FixAveCorrelate);
// clang-format on
#else

#ifndef LMP_FIX_AVE_CORRELATE_H
#define LMP_FIX_AVE_CORRELATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAveCorrelate : public Fix {
 public:
  FixAveCorrelate(class LAMMPS *, int, char **);
  ~FixAveCorrelate();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_array(int, int);

 private:
  int me, nvalues;
  int nrepeat, nfreq;
  bigint nvalid, nvalid_last;
  int *which, *argindex, *value2index;
  char **ids;
  FILE *fp;

  int type, ave, startstep, overwrite;
  double prefactor;
  long filepos;

  int firstindex;    // index in values ring of earliest time sample
  int lastindex;     // index in values ring of latest time sample
  int nsample;       // number of time samples in values ring

  int npair;    // number of correlation pairs to calculate
  int *count;
  double **values, **corr;

  int *save_count;    // saved values at Nfreq for output via compute_array()
  double **save_corr;

  void accumulate();
  bigint nextvalid();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open fix ave/correlate file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Compute ID for fix ave/correlate does not exist

Self-explanatory.

E: Fix ave/correlate compute does not calculate a scalar

Self-explanatory.

E: Fix ave/correlate compute does not calculate a vector

Self-explanatory.

E: Fix ave/correlate compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ID for fix ave/correlate does not exist

Self-explanatory.

E: Fix ave/correlate fix does not calculate a scalar

Self-explanatory.

E: Fix ave/correlate fix does not calculate a vector

Self-explanatory.

E: Fix ave/correlate fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix for fix ave/correlate not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/correlate
is requesting a value on a non-allowed timestep.

E: Variable name for fix ave/correlate does not exist

Self-explanatory.

E: Fix ave/correlate variable is not equal-style variable

Self-explanatory.

E: Fix ave/correlate variable is not vector-style variable

UNDOCUMENTED

E: Error writing file header

Something in the output to the file triggered an error.

E: Invalid timestep reset for fix ave/correlate

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Error writing out correlation data

Something in the output to the file triggered an error.

*/

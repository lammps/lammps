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

FixStyle(ave/correlate/long,FixAveCorrelateLong)

#else

#ifndef LMP_FIX_AVE_CORRELATE_LONG_H
#define LMP_FIX_AVE_CORRELATE_LONG_H

#include <cstdio>
#include "fix.h"

namespace LAMMPS_NS {

class FixAveCorrelateLong : public Fix {
 public:
  FixAveCorrelateLong(class LAMMPS *, int, char **);
  ~FixAveCorrelateLong();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

  void write_restart(FILE *);
  void restart(char *);
  double memory_usage();

  double *t; // Time steps for result arrays
  double **f; // Result arrays
  unsigned int npcorr;

 private:
  // NOT OPTIMAL: shift2 and accumulator2 only needed in cross-correlations
  double ***shift, *** shift2;
  double ***correlation;
  double **accumulator, **accumulator2;
  unsigned long int **ncorrelation;
  unsigned int *naccumulator;
  unsigned int *insertindex;

  int numcorrelators; // Recommended 20
  unsigned int p; // Points per correlator (recommended 16)
  unsigned int m; // Num points for average (recommended 2; p mod m = 0)
  unsigned int dmin; // Min distance between ponts for correlators k>0; dmin=p/m

  int length; // Length of result arrays
  int kmax; // Maximum correlator attained during simulation

  int me,nvalues;
  int nfreq;
  bigint nvalid,nvalid_last,last_accumulated_step;
  int *which,*argindex,*value2index;
  char **ids;
  FILE *fp;

  int type,startstep,overwrite;
  long filepos;

  int npair;           // number of correlation pairs to calculate
  double *values;

  void accumulate();
  void evaluate();
  bigint nextvalid();

  void add(const int i, const double w, const int k = 0);
  void add(const int i, const double wA, const double wB, const int k = 0);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open fix ave/correlate/long file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Compute ID for fix ave/correlate/long does not exist

Self-explanatory.

E: Fix ave/correlate/long compute does not calculate a scalar

Self-explanatory.

E: Fix ave/correlate/long compute does not calculate a vector

Self-explanatory.

E: Fix ave/correlate/long compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ID for fix ave/correlate/long does not exist

Self-explanatory.

E: Fix ave/correlate/long fix does not calculate a scalar

Self-explanatory.

E: Fix ave/correlate/long fix does not calculate a vector

Self-explanatory.

E: Fix ave/correlate/long fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix for fix ave/correlate/long not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/correlate/long
is requesting a value on a non-allowed timestep.

E: Variable name for fix ave/correlate/long does not exist

Self-explanatory.

E: Fix ave/correlate/long variable is not equal-style variable

Self-explanatory.

E: Invalid timestep reset for fix ave/correlate/long

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

*/

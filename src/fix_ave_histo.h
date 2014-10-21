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

FixStyle(ave/histo,FixAveHisto)

#else

#ifndef LMP_FIX_AVE_HISTO_H
#define LMP_FIX_AVE_HISTO_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAveHisto : public Fix {
 public:
  FixAveHisto(class LAMMPS *, int, char **);
  ~FixAveHisto();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_vector(int);
  double compute_array(int,int);
  void reset_timestep(bigint);

 private:
  int me,nvalues;
  int nrepeat,nfreq,irepeat;
  bigint nvalid;
  int *which,*argindex,*value2index;
  char **ids;
  FILE *fp;
  double lo,hi,binsize,bininv;
  int kind,beyond,overwrite;
  long filepos;

  double stats[4],stats_total[4],stats_all[4];
  double **stats_list;

  int nbins;
  double *bin,*bin_total,*bin_all;
  double **bin_list;
  double *coord;

  double *vector;
  int maxatom;

  int ave,nwindow,startstep,mode;
  char *title1,*title2,*title3;
  int iwindow,window_limit;

  void bin_one(double);
  void bin_vector(int, double *, int);
  void bin_atoms(double *, int);
  void options(int, char **);
  void allocate_values(int);
  bigint nextvalid();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for fix ave/histo does not exist

Self-explanatory.

E: Fix ID for fix ave/histo does not exist

Self-explanatory.

E: Fix ave/histo input is invalid compute

Self-explanatory.

E: Fix ave/histo input is invalid fix

Self-explanatory.

E: Fix ave/histo input is invalid variable

Self-explanatory.

E: Fix ave/histo inputs are not all global, peratom, or local

All inputs in a single fix ave/histo command must be of the
same style.

E: Fix ave/histo cannot input per-atom values in scalar mode

Self-explanatory.

E: Fix ave/histo cannot input local values in scalar mode

Self-explanatory.

E: Fix ave/histo compute does not calculate a global scalar

Self-explanatory.

E: Fix ave/histo compute does not calculate a global vector

Self-explanatory.

E: Fix ave/histo compute vector is accessed out-of-range

Self-explanatory.

E: Fix ave/histo compute does not calculate a global array

Self-explanatory.

E: Fix ave/histo compute array is accessed out-of-range

Self-explanatory.

E: Fix ave/histo compute does not calculate per-atom values

Self-explanatory.

E: Fix ave/histo compute does not calculate a per-atom vector

Self-explanatory.

E: Fix ave/histo compute does not calculate a per-atom array

Self-explanatory.

E: Fix ave/histo compute does not calculate local values

Self-explanatory.

E: Fix ave/histo compute does not calculate a local vector

Self-explanatory.

E: Fix ave/histo compute does not calculate a local array

Self-explanatory.

E: Fix ave/histo fix does not calculate a global scalar

Self-explanatory.

E: Fix ave/histo fix does not calculate a global vector

Self-explanatory.

E: Fix ave/histo fix vector is accessed out-of-range

Self-explanatory.

E: Fix for fix ave/histo not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/histo is
requesting a value on a non-allowed timestep.

E: Fix ave/histo fix does not calculate a global array

Self-explanatory.

E: Fix ave/histo fix array is accessed out-of-range

Self-explanatory.

E: Fix ave/histo fix does not calculate per-atom values

Self-explanatory.

E: Fix ave/histo fix does not calculate a per-atom vector

Self-explanatory.

E: Fix ave/histo fix does not calculate a per-atom array

Self-explanatory.

E: Fix ave/histo fix does not calculate local values

Self-explanatory.

E: Fix ave/histo fix does not calculate a local vector

Self-explanatory.

E: Fix ave/histo fix does not calculate a local array

Self-explanatory.

E: Variable name for fix ave/histo does not exist

Self-explanatory.

E: Cannot open fix ave/histo file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Fix ave/histo missed timestep

You cannot reset the timestep to a value beyond where the fix
expects to next perform averaging.

*/

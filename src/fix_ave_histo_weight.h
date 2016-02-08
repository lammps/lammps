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

FixStyle(ave/histo/weight,FixAveHistoWeight)

#else

#ifndef LMP_FIX_AVE_HISTO_WEIGHT_H
#define LMP_FIX_AVE_HISTO_WEIGHT_H

#include <stdio.h>
#include "fix_ave_histo.h"

namespace LAMMPS_NS {

class FixAveHistoWeight : public FixAveHisto {
 public:
  FixAveHistoWeight(class LAMMPS *, int, char **);
  ~FixAveHistoWeight() {}
  void end_of_step();

 private:
  void bin_one_weights(double, double);
  void bin_vector_weights(int, double *, int, double *, int);
  void bin_atoms_weights(double *, int, double *, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix ave/histo/weight value and weight vector lengths do not match

UNDOCUMENTED

E: Invalid timestep reset for fix ave/histo

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Error writing out histogram data

UNDOCUMENTED

U: Compute ID for fix ave/histo does not exist

Self-explanatory.

U: Fix ID for fix ave/histo does not exist

Self-explanatory.

U: Fix ave/histo input is invalid compute

Self-explanatory.

U: Fix ave/histo input is invalid fix

Self-explanatory.

U: Fix ave/histo input is invalid variable

Self-explanatory.

U: Fix ave/histo inputs are not all global, peratom, or local

All inputs in a single fix ave/histo command must be of the
same style.

U: Fix ave/histo cannot input per-atom values in scalar mode

Self-explanatory.

U: Fix ave/histo cannot input local values in scalar mode

Self-explanatory.

U: Fix ave/histo compute does not calculate a global scalar

Self-explanatory.

U: Fix ave/histo compute does not calculate a global vector

Self-explanatory.

U: Fix ave/histo compute vector is accessed out-of-range

Self-explanatory.

U: Fix ave/histo compute does not calculate a global array

Self-explanatory.

U: Fix ave/histo compute array is accessed out-of-range

Self-explanatory.

U: Fix ave/histo compute does not calculate per-atom values

Self-explanatory.

U: Fix ave/histo compute does not calculate a per-atom vector

Self-explanatory.

U: Fix ave/histo compute does not calculate a per-atom array

Self-explanatory.

U: Fix ave/histo compute does not calculate local values

Self-explanatory.

U: Fix ave/histo compute does not calculate a local vector

Self-explanatory.

U: Fix ave/histo compute does not calculate a local array

Self-explanatory.

U: Fix ave/histo fix does not calculate a global scalar

Self-explanatory.

U: Fix ave/histo fix does not calculate a global vector

Self-explanatory.

U: Fix ave/histo fix vector is accessed out-of-range

Self-explanatory.

U: Fix for fix ave/histo not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/histo is
requesting a value on a non-allowed timestep.

U: Fix ave/histo fix does not calculate a global array

Self-explanatory.

U: Fix ave/histo fix array is accessed out-of-range

Self-explanatory.

U: Fix ave/histo fix does not calculate per-atom values

Self-explanatory.

U: Fix ave/histo fix does not calculate a per-atom vector

Self-explanatory.

U: Fix ave/histo fix does not calculate a per-atom array

Self-explanatory.

U: Fix ave/histo fix does not calculate local values

Self-explanatory.

U: Fix ave/histo fix does not calculate a local vector

Self-explanatory.

U: Fix ave/histo fix does not calculate a local array

Self-explanatory.

U: Variable name for fix ave/histo does not exist

Self-explanatory.

U: Cannot open fix ave/histo file %s

The specified file cannot be opened.  Check that the path and name are
correct.

*/

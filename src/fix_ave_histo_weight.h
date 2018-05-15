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

#include <cstdio>
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

Self-explanatory.

E: Invalid timestep reset for fix ave/histo

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Fix ave/histo/weight option not yet supported

UNDOCUMENTED

E: Error writing out histogram data

Something in the output to the file triggered an error.

*/

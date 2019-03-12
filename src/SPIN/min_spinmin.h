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

#ifdef MINIMIZE_CLASS

MinimizeStyle(spinmin,MinSpinMin)

#else

#ifndef LMP_MIN_SPINMIN_H
#define LMP_MIN_SPINMIN_H

#include "min.h"

namespace LAMMPS_NS {

class MinSpinMin : public Min {
 public:
  MinSpinMin(class LAMMPS *);
  ~MinSpinMin() {}
  void init();
  void setup_style();
  void reset_vectors();
  int iterate(int);
  double evaluate_dt();
  void advance_spins(double);
  double fmnorm_sqr();

 private:

  // global and spin timesteps
  
  double dt;
  double dts;

  double *spvec;               // variables for atomic dof, as 1d vector
  double *fmvec;               // variables for atomic dof, as 1d vector

  bigint last_negative;
};

}

#endif
#endif

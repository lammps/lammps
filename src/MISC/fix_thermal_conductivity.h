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

FixStyle(thermal/conductivity,FixThermalConductivity)

#else

#ifndef LMP_FIX_THERMAL_CONDUCTIVITY_H
#define LMP_FIX_THERMAL_CONDUCTIVITY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixThermalConductivity : public Fix {
 public:
  FixThermalConductivity(class LAMMPS *, int, char **);
  ~FixThermalConductivity();
  int setmask();
  void init();
  void end_of_step();
  double compute_scalar();

 private:
  int me;
  int edim,nbin,periodicity;
  int nswap;
  double prd,boxlo,boxhi;
  double slablo_lo,slablo_hi,slabhi_lo,slabhi_hi;
  double e_exchange;

  int nlo,nhi;
  int *index_lo,*index_hi;
  double *ke_lo,*ke_hi;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix thermal/conductivity swap value must be positive

Self-explanatory.

W: Fix thermal/conductivity comes before fix ave/spatial

The order of these 2 fixes in your input script is such that fix
thermal/conductivity comes first.  If you are using fix ave/spatial to
measure the temperature profile induced by fix viscosity, then this
may cause a glitch in the profile since you are averaging immediately
after swaps have occurred.  Flipping the order of the 2 fixes
typically helps.

*/

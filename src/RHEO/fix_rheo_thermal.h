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

FixStyle(rheo/thermal,FixRHEOThermal)

#else

#ifndef LMP_FIX_RHEO_THERMAL_H
#define LMP_FIX_RHEO_THERMAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRHEOThermal : public Fix {
 public:
  FixRHEOThermal(class LAMMPS *, int, char **);
  ~FixRHEOThermal();
  int setmask();
  void init();
  void setup_pre_force(int);
  void initial_integrate();
  void post_integrate();
  void pre_force(int);
  void end_of_step();
  void reset_dt();

 private:
  double *cv_type, cv;
  double *Tc_type, Tc;
  double *kappa_type, kappa;
  double *alpha_type, alpha;
  double dtf, dtv;

  double calc_kappa(int);
  double calc_cv(int);
  double calc_Tc(int);
  double calc_alpha(int);

  int Tc_style;
  int cv_style;
  int alpha_style;
  int conductivity_style;

  class FixRHEO *fix_rheo;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/

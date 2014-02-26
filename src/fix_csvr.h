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

FixStyle(csvr,FixCSVR)

#else

#ifndef LMP_FIX_CSVR_H
#define LMP_FIX_CSVR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCSVR : public Fix {
 public:
  FixCSVR(class LAMMPS *, int, char **);
  virtual ~FixCSVR() {}

  int setmask();
  int modify_param(int, char **);
  void reset_target(double);

  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void initial_integrate_respa(int, int, int);
  virtual void final_integrate_respa(int, int);
  virtual void reset_dt();
  
  double compute_scalar();
  virtual void *extract(const char *, int &);

 protected:
  class RanMars *random;

  double t_start,t_stop,t_period,t_target;
  double dtv,dtf;
  double *step_respa;
  int mass_require;
  int which;

  double energy;
  int tstyle,tvar;
  char *tstr;

  char *id_temp;
  class Compute *temperature;
  int tflag;
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

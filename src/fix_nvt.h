/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_NVT_H
#define FIX_NVT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNVT : public Fix {
 public:
  FixNVT(class LAMMPS *, int, char **);
  virtual ~FixNVT();
  int setmask();
  virtual void init();
  void setup(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int);
  double compute_scalar();
  void write_restart(FILE *);
  void restart(char *);
  int modify_param(int, char **);
  void reset_target(double);
  void reset_dt();

 protected:
  int which,chain;
  double t_start,t_stop;
  double t_current,t_target;
  double t_freq,drag,drag_factor;
  double f_eta,eta_dot,eta,eta2_dot,eta2,factor,factor2;
  double dtv,dtf,dthalf,dt4,dt8;

  int nlevels_respa;
  double *step_respa;

  char *id_temp;
  class Compute *temperature;
  int tflag;
};

}

#endif

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_NVT_H
#define FIX_NVT_H

#include "fix.h"
class Temperature;

class FixNVT : public Fix {
 public:
  FixNVT(int, char **);
  ~FixNVT() {}
  int setmask();
  void init();
  void setup();
  void initial_integrate();
  void final_integrate();
  void initial_integrate_respa(int,int);
  void final_integrate_respa(int);
  void write_restart(FILE *);
  void restart(char *);
  int modify_param(int, char **);
  void reset_target(double);
  int thermo_fields(int, int *, char **);
  int thermo_compute(double *);

 private:
  double t_start,t_stop;
  double t_current,t_target;
  double t_freq,drag,drag_factor;
  double f_eta,eta_dot,eta,factor;
  double dtv,dtf,dthalf;
  Temperature *temperature;

  int nlevels_respa;
  double *step_respa;
};

#endif

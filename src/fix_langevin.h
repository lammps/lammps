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

#ifndef FIX_LANGEVIN_H
#define FIX_LANGEVIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevin : public Fix {
 public:
  FixLangevin(class LAMMPS *, int, char **);
  ~FixLangevin();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void end_of_step();
  void reset_target(double);
  void reset_dt();
  int modify_param(int, char **);
  double compute_scalar();
  double memory_usage();

 private:
  int which,tally;
  double t_start,t_stop,t_period;
  double *gfactor1,*gfactor2,*ratio;
  double energy,energy_onestep;

  int nmax;
  double **flangevin;

  char *id_temp;
  class Compute *temperature;

  int nlevels_respa;
  class RanMars *random;

  void post_force_no_tally();
  void post_force_tally();
};

}

#endif

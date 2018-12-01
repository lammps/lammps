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

FixStyle(ffl,FixFFL)

#else

#ifndef LMP_FIX_FFL_H
#define LMP_FIX_FFL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFFL : public Fix {
 public:
  FixFFL(class LAMMPS *, int, char **);
  virtual ~FixFFL();
  int setmask();
  void init();
  void setup(int);
  void ffl_integrate();
  void initial_integrate_respa(int vflag, int ilevel, int iloop);
  void final_integrate_respa(int ilevel, int iloop);
  void initial_integrate(int vflag);
  void final_integrate();
  double compute_scalar();
  void reset_target(double);
  virtual void reset_dt();
  double memory_usage();
  void grow_arrays(int);

  virtual void *extract(const char *, int &);

  void init_ffl();
 protected:
  double *ffl_tmp1, *ffl_tmp2;
  double t_start, t_stop, t_target;
  double dtv, dtf, c1, c2, gamma;
  char flip_type[10];

  int doffl, ffl_every, ffl_step, flip_int;
  class RanMars *random;
  double *sqrt_m;
  double *step_respa;
  double energy;
  int nlevels_respa;

  double **vaux;
};

}

#endif
#endif

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

#ifndef FIX_WALL_LJ126_H
#define FIX_WALL_LJ126_H

#include "fix.h"

class FixWallLJ126 : public Fix {
 public:
  FixWallLJ126(int, char **);
  ~FixWallLJ126() {}
  int setmask();
  void init();
  void setup();
  void min_setup();
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);

  int thermo_fields(int, int *, char **);
  int thermo_compute(double *);

 private:
  int dim,side,thermo_flag,eflag_on;
  double coord,epsilon,sigma,cutoff;
  double coeff1,coeff2,coeff3,coeff4,offset;
  double eng,etotal;
  int nlevels_respa;
};

#endif


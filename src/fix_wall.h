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

#ifndef LMP_FIX_WALL_H
#define LMP_FIX_WALL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWall : public Fix {
 public:
  FixWall(class LAMMPS *, int, char **);
  virtual ~FixWall() {}
  int setmask();
  virtual void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

  virtual void precompute(int) = 0;
  virtual void wall_particle(int, double) = 0;

 protected:
  int wallflag[6];
  double coord0[6],epsilon[6],sigma[6],cutoff[6];
  int velflag,wigflag;
  double vel[6],amplitude[6];
  double period,omega;
  int eflag;
  double ewall[7],ewall_all[7];
  int nlevels_respa;
  double dt;
  int time_origin;
};

}

#endif

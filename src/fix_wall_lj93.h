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

#ifndef FIX_WALL_LJ93_H
#define FIX_WALL_LJ93_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallLJ93 : public Fix {
 public:
  FixWallLJ93(class LAMMPS *, int, char **);
  ~FixWallLJ93() {}
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  int dim,side,thermo_flag,eflag_enable;
  double coord,epsilon,sigma,cutoff;
  double coeff1,coeff2,coeff3,coeff4,offset;
  double wall[4],wall_all[4];
  double vel,coord0;
  int wall_flag;
  int nlevels_respa;
};

}

#endif

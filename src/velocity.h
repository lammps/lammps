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

#ifndef VELOCITY_H
#define VELOCITY_H

#include "pointers.h"

namespace LAMMPS_NS {

class Velocity : protected Pointers {
 public:
  Velocity(class LAMMPS *);
  void command(int, char **);
  void init_external(char *);
  void options(int, char **);
  void create(double, int);

 private:
  int igroup,groupbit;
  int style;
  int dist_flag,sum_flag,momentum_flag,rotation_flag,loop_flag,scale_flag;
  double xscale,yscale,zscale;
  class Compute *temperature;

  void set(int, char **);
  void scale(int, char **);
  void ramp(int, char **);
  void zero(int, char **);

  void rescale(double, double);
  void zero_momentum();
  void zero_rotation();
};

}

#endif

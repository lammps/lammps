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

#ifndef VELOCITY_H
#define VELOCITY_H

#include "lammps.h"

class Temperature;
class RanPark;

class Velocity : public LAMMPS {
 public:
  Velocity() {}
  ~Velocity() {}
  void command(int, char **);

 private:
  int igroup,groupbit;
  int style,tempwhich;
  int dist_flag,sum_flag,momentum_flag,rotation_flag,loop_flag,scale_flag;
  double xscale,yscale,zscale;
  Temperature *temperature;

  void create(int, char **);
  void set(int, char **);
  void scale(int, char **);
  void ramp(int, char **);
  void zero(int, char **);

  void rescale(double, double);
  void zero_momentum();
  void zero_rotation();

  void options(int, char **);
  void triple(double, double, double, double *, double *, double *,
	      int, RanPark *);
};

#endif

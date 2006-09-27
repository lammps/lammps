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

#ifndef ANGLE_H
#define ANGLE_H

#include "stdio.h"
#include "lammps.h"

class Angle : public LAMMPS {
 public:
  int allocated;
  int *setflag;
  double energy;
  double virial[6];
  double PI;

  Angle();
  virtual ~Angle() {}
  virtual void init();
  virtual void compute(int, int) = 0;
  virtual void settings(int, char **) {}
  virtual void coeff(int, int, char **) = 0;
  virtual double equilibrium_angle(int) = 0;
  virtual void write_restart(FILE *) = 0;
  virtual void read_restart(FILE *) = 0;
  virtual double single(int, int, int, int, double) {return 0.0;}
  virtual int memory_usage() {return 0;}
};

#endif

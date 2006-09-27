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

#ifndef ANGLE_CLASS2_H
#define ANGLE_CLASS2_H

#include "stdio.h"
#include "angle.h"

class AngleClass2 : public Angle {
 public:
  AngleClass2() {}
  ~AngleClass2();
  void compute(int, int);
  void coeff(int, int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
  double *theta0,*k2,*k3,*k4;
  double *bb_k,*bb_r1,*bb_r2;
  double *ba_k1,*ba_k2,*ba_r1,*ba_r2;
  int *setflag_a,*setflag_bb,*setflag_ba;

  void allocate();
};

#endif

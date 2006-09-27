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

#ifndef REGION_SPHERE_H
#define REGION_SPHERE_H

#include "region.h"

class RegSphere : public Region {
 public:
  RegSphere(int, char **);
  ~RegSphere() {}
  int match(double, double, double);

 private:
  double xc,yc,zc;
  double radius;
};

#endif

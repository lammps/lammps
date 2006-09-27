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

#ifndef REGION_CYLINDER_H
#define REGION_CYLINDER_H

#include "region.h"

class RegCylinder : public Region {
  friend class FixInsert;

 public:
  RegCylinder(int, char **);
  ~RegCylinder() {}
  int match(double, double, double);

 private:
  char axis;
  double c1,c2;
  double radius;
  double lo,hi;
};

#endif

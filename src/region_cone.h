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

#ifdef REGION_CLASS

RegionStyle(cone,RegCone)

#else

#ifndef LMP_REGION_CONE_H
#define LMP_REGION_CONE_H

#include "region.h"

namespace LAMMPS_NS {

class RegCone : public Region {
 public:
  RegCone(class LAMMPS *, int, char **);
  ~RegCone();
  int inside(double, double, double);
  int surface_interior(double *, double);
  int surface_exterior(double *, double);

 private:
  char axis;
  double c1,c2;
  double radiuslo,radiushi;
  double lo,hi;
  double maxradius;

  void point_on_line_segment(double *, double *, double *, double *);
  double closest(double *, double *, double *, double);

  void subtract(double *, double *, double *);
  double dotproduct(double *, double *);
};

}

#endif
#endif

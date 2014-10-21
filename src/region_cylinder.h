/* -*- c++ -*- ----------------------------------------------------------
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

RegionStyle(cylinder,RegCylinder)

#else

#ifndef LMP_REGION_CYLINDER_H
#define LMP_REGION_CYLINDER_H

#include "region.h"

namespace LAMMPS_NS {

class RegCylinder : public Region {
  friend class FixPour;

 public:
  RegCylinder(class LAMMPS *, int, char **);
  ~RegCylinder();
  void init();
  int inside(double, double, double);
  int surface_interior(double *, double);
  int surface_exterior(double *, double);
  void shape_update();

 private:
  char axis;
  double c1,c2;
  double radius;
  double lo,hi;
  int rstyle,rvar;
  char *rstr;

  void variable_check();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use region INF or EDGE when box does not exist

Regions that extend to the box boundaries can only be used after the
create_box command has been used.

E: Variable evaluation in region gave bad value

Variable returned a radius < 0.0.

E: Variable name for region cylinder does not exist

Self-explanatory.

E: Variable for region cylinder is invalid style

Only equal-style varaibles are allowed.

*/

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

RegionStyle(sphere,RegSphere)

#else

#ifndef LMP_REGION_SPHERE_H
#define LMP_REGION_SPHERE_H

#include "region.h"

namespace LAMMPS_NS {

class RegSphere : public Region {
 public:
  RegSphere(class LAMMPS *, int, char **);
  ~RegSphere();
  void init();
  int inside(double, double, double);
  int surface_interior(double *, double);
  int surface_exterior(double *, double);
  void shape_update();
  void set_velocity_shape();
  void velocity_contact_shape(double *, double *);


 private:
  double xc,yc,zc;
  double radius;
  int xstyle,xvar;
  int ystyle,yvar;
  int zstyle,zvar;
  int rstyle,rvar;
  char *xstr,*ystr,*zstr,*rstr;

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

E: Variable evaluation in region gave bad value

Variable returned a radius < 0.0.

E: Variable name for region sphere does not exist

Self-explanatory.

E: Variable for region sphere is invalid style

Only equal-style variables are allowed.

*/

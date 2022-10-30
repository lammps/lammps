/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS
// clang-format off
RegionStyle(sphere,RegSphere);
// clang-format on
#else

#ifndef LMP_REGION_SPHERE_H
#define LMP_REGION_SPHERE_H

#include "region.h"

namespace LAMMPS_NS {

class RegSphere : public Region {
 public:
  RegSphere(class LAMMPS *, int, char **);
  ~RegSphere() override;
  void init() override;
  int inside(double, double, double) override;
  int surface_interior(double *, double) override;
  int surface_exterior(double *, double) override;
  void shape_update() override;
  void set_velocity_shape() override;
  void velocity_contact_shape(double *, double *) override;

 private:
  double xc, yc, zc;
  double radius;
  int xstyle, xvar;
  int ystyle, yvar;
  int zstyle, zvar;
  int rstyle, rvar;
  char *xstr, *ystr, *zstr, *rstr;

  void variable_check();
};

}    // namespace LAMMPS_NS

#endif
#endif

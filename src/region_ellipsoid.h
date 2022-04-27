/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS
// clang-format off
RegionStyle(ellipsoid,RegEllipsoid);
// clang-format on
#else

#ifndef LMP_REGION_ELLIPSOID_H
#define LMP_REGION_ELLIPSOID_H

#include "region.h"

namespace LAMMPS_NS {

class RegEllipsoid : public Region {
 public:
  RegEllipsoid(class LAMMPS *, int, char **);
  ~RegEllipsoid() override;
  void init() override;
  int inside(double, double, double) override;
  int surface_interior(double *, double) override;
  int surface_exterior(double *, double) override;
  void shape_update() override;

 private:
  double xc, yc, zc;
  double a, b, c;
  int xstyle, xvar;
  int ystyle, yvar;
  int zstyle, zvar;
  int astyle, avar;
  int bstyle, bvar;
  int cstyle, cvar;
  char *xstr, *ystr, *zstr;
  char *astr, *bstr, *cstr;

  void variable_check();

  double GetRoot2D(double r0, double z0, double z1, double g);
  double GetRoot3D(double r0, double r1, double z0, double z1, double z2, double g);
  double DistancePointEllipse(double e0, double e1, double y0, double y1, double &x0, double &x1);
  double DistancePointEllipsoid(double e0, double e1, double e2, double y0, double y1, double y2,
                                double &x0, double &x1, double &x2);
};

}    // namespace LAMMPS_NS

#endif
#endif

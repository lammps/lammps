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
};

}    // namespace LAMMPS_NS

#endif
#endif

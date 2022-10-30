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
RegionStyle(cone,RegCone);
// clang-format on
#else

#ifndef LMP_REGION_CONE_H
#define LMP_REGION_CONE_H

#include "region.h"

namespace LAMMPS_NS {

class RegCone : public Region {
 public:
  RegCone(class LAMMPS *, int, char **);
  ~RegCone() override;
  int inside(double, double, double) override;
  int surface_interior(double *, double) override;
  int surface_exterior(double *, double) override;

 private:
  char axis;
  double c1, c2;
  double radiuslo, radiushi;
  double lo, hi;
  double maxradius;

  double closest(double *, double *, double *, double);
};

}    // namespace LAMMPS_NS

#endif
#endif

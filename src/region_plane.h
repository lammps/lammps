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
RegionStyle(plane,RegPlane);
// clang-format on
#else

#ifndef LMP_REGION_PLANE_H
#define LMP_REGION_PLANE_H

#include "region.h"

namespace LAMMPS_NS {

class RegPlane : public Region {
 public:
  RegPlane(class LAMMPS *, int, char **);
  ~RegPlane() override;
  int inside(double, double, double) override;
  int surface_interior(double *, double) override;
  int surface_exterior(double *, double) override;

 private:
  double xp, yp, zp;
  double normal[3];
};

}    // namespace LAMMPS_NS

#endif
#endif

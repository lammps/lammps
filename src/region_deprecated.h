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
// list all deprecated and removed region styles here
RegionStyle(DEPRECATED,RegionDeprecated);
// clang-format on
#else

#ifndef LMP_REGION_DEPRECATED_H
#define LMP_REGION_DEPRECATED_H

#include "region.h"

namespace LAMMPS_NS {

class RegionDeprecated : public Region {
 public:
  RegionDeprecated(class LAMMPS *, int, char **);

  void init() override {}
  int inside(double, double, double) override { return 0; }
  int surface_interior(double *, double) override { return 0; }
  int surface_exterior(double *, double) override { return 0; }
};

}    // namespace LAMMPS_NS

#endif
#endif

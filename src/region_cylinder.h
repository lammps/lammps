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
RegionStyle(cylinder,RegCylinder);
// clang-format on
#else

#ifndef LMP_REGION_CYLINDER_H
#define LMP_REGION_CYLINDER_H

#include "region.h"

namespace LAMMPS_NS {

class RegCylinder : public Region {
  friend class FixPour;

 public:
  RegCylinder(class LAMMPS *, int, char **);
  ~RegCylinder() override;
  void init() override;
  int inside(double, double, double) override;
  int surface_interior(double *, double) override;
  int surface_exterior(double *, double) override;
  void shape_update() override;
  void set_velocity_shape() override;
  void velocity_contact_shape(double *, double *) override;

 private:
  char axis;
  double c1, c2;
  double radius;
  double lo, hi;
  int c1style, c1var;
  int c2style, c2var;
  int rstyle, rvar;
  char *c1str, *c2str, *rstr;

  void variable_check();
};

}    // namespace LAMMPS_NS

#endif
#endif

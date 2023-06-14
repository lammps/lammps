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
RegionStyle(block,RegBlock);
// clang-format on
#else

#ifndef LMP_REGION_BLOCK_H
#define LMP_REGION_BLOCK_H

#include "region.h"

namespace LAMMPS_NS {

class RegBlock : public Region {
  friend class FixPour;

 public:
  RegBlock(class LAMMPS *, int, char **);
  ~RegBlock() override;
  void init() override;
  int inside(double, double, double) override;
  int surface_interior(double *, double) override;
  int surface_exterior(double *, double) override;
  void shape_update() override;

 protected:
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double corners[6][4][3];
  double face[6][3];
  int xlostyle, xlovar, xhistyle, xhivar;
  int ylostyle, ylovar, yhistyle, yhivar;
  int zlostyle, zlovar, zhistyle, zhivar;
  char *xlostr, *ylostr, *zlostr;
  char *xhistr, *yhistr, *zhistr;

  double find_closest_point(int, double *, double &, double &, double &);
  int inside_face(double *, int);
  void variable_check();
};

}    // namespace LAMMPS_NS

#endif
#endif

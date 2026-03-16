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
RegionStyle(grid,RegGrid);
// clang-format on
#else

#ifndef LMP_REGION_GRID_H
#define LMP_REGION_GRID_H

#include "region.h"

namespace LAMMPS_NS {

class RegGrid : public Region {
 public:
  RegGrid(class LAMMPS *, int, char **);
  ~RegGrid() override;
  void init() override;
  void shape_update() override;
  int inside(double, double, double) override;
  int surface_interior(double *, double) override;
  int surface_exterior(double *, double) override;

 protected:
  enum { COMPUTE_SOURCE, FIX_SOURCE };
  enum { OP_GT, OP_GE, OP_LT, OP_LE, OP_EQ, OP_NE };

  char *gridref;
  int source_type;
  char *source_id;
  char *grid_name;
  char *data_name;
  int igrid, idata;
  int gridindex;
  int compare_op;
  double threshold;

  class Grid3d *grid3d;
  int nx, ny, nz, ncol;
  int nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out;
  double ***vec3d, ****array3d;

  void resolve_grid_reference();
  void update_bbox();
  int evaluate(double value);
  bool cell_inside(int, int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif

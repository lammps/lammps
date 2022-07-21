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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(property/grid,ComputePropertyGrid);
// clang-format on
#else

#ifndef LMP_COMPUTE_PROPERTY_GRID_H
#define LMP_COMPUTE_PROPERTY_GRID_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePropertyGrid : public Compute {
 public:
  ComputePropertyGrid(class LAMMPS *, int, char **);
  ~ComputePropertyGrid() override;
  void init() override {}
  void compute_pergrid() override;

  int get_grid_by_name(char *, int &);
  void *get_grid_by_index(int);
  int get_griddata_by_name(int, char *, int &);
  void *get_griddata_by_index(int);

  double memory_usage() override;

 private:
  int nx,ny,nz;
  int nvalues;
  int dimension;

  class Grid3d *gc;
  int ngc_buf1, ngc_buf2;
  double *gc_buf1, *gc_buf2;

  double ****data;

  typedef void (ComputePropertyGrid::*FnPtrPack)(int);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_id(int);

  void pack_x(int);
  void pack_y(int);
  void pack_z(int);

  void pack_xs(int);
  void pack_ys(int);
  void pack_zs(int);

  void pack_xc(int);
  void pack_yc(int);
  void pack_zc(int);

  void pack_xsc(int);
  void pack_ysc(int);
  void pack_zsc(int);
};

}    // namespace LAMMPS_NS

#endif
#endif

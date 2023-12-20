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
  void init() override;
  void compute_pergrid() override;

  void reset_grid() override;

  int get_grid_by_name(const std::string &, int &) override;
  void *get_grid_by_index(int) override;
  int get_griddata_by_name(int, const std::string &, int &) override;
  void *get_griddata_by_index(int) override;

  double memory_usage() override;

 private:
  int nxgrid,nygrid,nzgrid;
  int nvalues;
  int dimension;
  int triclinic;

  class Grid2d *grid2d;
  class Grid3d *grid3d;

  int nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in;
  int nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out;
  int ngridout;

  double **vec2d,***vec3d;
  double ***array2d,****array3d;

  // local methods

  void allocate_grid();
  void deallocate_grid();

  typedef void (ComputePropertyGrid::*FnPtrPack)(int);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_id(int);
  void pack_proc(int);
  template <int IDIM> void pack_indices(int);
  template <int POS, int MODE, int IDIM> void pack_coords(int);
};

}    // namespace LAMMPS_NS

#endif
#endif

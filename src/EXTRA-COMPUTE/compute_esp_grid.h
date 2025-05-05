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
ComputeStyle(esp/grid, ComputeESPGrid)
// clang-format on
#else

#ifndef LMP_COMPUTE_ESP_GRID_H
#define LMP_COMPUTE_ESP_GRID_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeESPGrid : public Compute {
 public:
  ComputeESPGrid(class LAMMPS *, int, char **);
  ~ComputeESPGrid() override;

  void init() override;
  double compute_scalar() override;
  void compute_vector() override;
  void compute_pergrid() override;
  void reset_grid() override;

  int get_grid_by_name(const std::string &, int &) override;
  void *get_grid_by_index(int) override;
  int get_griddata_by_name(int, const std::string &, int &) override;
  void *get_griddata_by_index(int) override;

 private:
  double spacing;
  int nx, ny, nz;

  int ixlo, ixhi, iylo, iyhi, izlo, izhi;
  int oxlo, oxhi, oylo, oyhi, ozlo, ozhi;

  double xlo, ylo, zlo;

  class Grid3d *esp_grid;
  double ***esp, ***reference, ***weight;
  double *bcut_acks2;
  int reaxflag;

  void allocate_grid();
  void deallocate_grid();

  inline double compute_weight(double r, double rcut) const;

};

}  // namespace LAMMPS_NS

#endif
#endif


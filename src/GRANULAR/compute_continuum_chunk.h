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
ComputeStyle(continuum/chunk,ComputeContinuumChunk);
// clang-format on
#else

#ifndef LMP_COMPUTE_CONTINUUM_CHUNK_H
#define LMP_COMPUTE_CONTINUUM_CHUNK_H

#include "compute_chunk.h"

#include <string>
#include <vector>

namespace LAMMPS_NS {

class ComputeContinuumChunk : public ComputeChunk {
 public:
  ComputeContinuumChunk(class LAMMPS *, int, char **);
  ~ComputeContinuumChunk() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  double memory_usage() override;

 private:
  std::vector<std::pair<int, int>> values;
  std::vector<std::string> labels;

  struct StencilOffset {
    int dn[3];     // index shifts along chunk axes 0, 1, 2
    double dx[3];  // displacement shifts along spatial axes x, y, z
  };
  std::vector<StencilOffset> stencil;

  int dim, bin_dim, pstyle, calculate_pair, calculate_2_loops;
  int boundary_group_flag, boundary_groupbit;
  int index_density, index_momentum[3], index_velocity[3], index_vgrad[3][3];
  double w_cut, w_cut_sq, w_sd, w_sd_sq, w_scale, w_offset;

  int nvalues, nskip, radius_required;
  int boundaryflag;

  class NeighList *list;

  double *delta;
  int ncoord, reducedflag;
  int *nlayers, *chunk_dim;

  double **values_local, **values_global;
  double *density_local, *density_global;
  double **momentum_local, **momentum_global;

  void allocate() override;
  inline double calc_w(double) const;
  inline double calc_w_int(double *, double *) const;
  void add_tensor_component(char *, int);
  void add_vector_component(char *, int);
  int shifted_bin(int, int *) const;
  void build_stencil();
  std::string get_thermo_colname(int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.
   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
   the U.S. Government retains certain rights in this software.

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// Contributing author, Richard Meng, Queen's University at Kingston, 22.11.24,
// contact@richardzjm.com

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mtp,PairMTP);
// clang-format on
#else

#ifndef LMP_PAIR_MTP_H
#define LMP_PAIR_MTP_H

#include "mtp_radial_basis.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairMTP : public Pair {
 public:
  PairMTP(class LAMMPS *);
  ~PairMTP() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

 protected:
  // Atoms per pass. Used only for Kokkos versions. Processed here for CPU parsing.
  static constexpr int DEFAULT_CHUNKSIZE = 65536;
  int input_chunk_size;

  virtual int settings_keyword(int narg, char **arg, int iarg);
  virtual void allocate();
  virtual void read_file(FILE *);
  void prepare_map(int narg, char **arg);
  void prepare_angular();

  std::string potential_name;
  std::string potential_tag;
  int species_count;
  double scaling;

  // Radial basis: 1 => RBChebyshev.
  int radial_basis_type_index;
  RadialMTPBasis *radial_basis;
  double *radial_basis_coeffs;
  int radial_func_count;
  int radial_basis_size;
  int radial_coeff_count;
  int radial_coeff_count_per_pair;

  // One cutoff for all species combinations.
  double min_cutoff;
  double max_cutoff;
  double max_cutoff_sq;

  double *linear_coeffs;
  double *species_coeffs;
  int alpha_moment_count, alpha_index_basic_count, alpha_index_times_count, alpha_scalar_count,
      max_alpha_index_basic;
  int **alpha_index_basic;
  int **alpha_index_times;
  int force_index_times_count;
  int *alpha_moment_mapping;

  // Shared angular monomials and per-neighbor scratch.
  int angular_count;
  int *basic_to_angular;
  int *basic_by_mu;
  int *angular_by_mu;
  int *mu_offsets;
  int *angular_parent;
  int *angular_axis;
  double *angular_vals;
  double *angular_ders;
  double *basic_ders_by_mu;

  // Graph traversal, forwards and backwards pass
  double *moment_tensor_vals;
  double *nbh_energy_ders_wrt_moments;

  // Cache values between forwards and backwards pass
  int cache_size;
  int *cached_j;
  double **neighbor_cache;
};

}    // namespace LAMMPS_NS

#endif
#endif

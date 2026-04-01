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

//
// Contributing author, Richard Meng, Queen's University at Kingston, 22.11.24, contact@richardzjm.com
//

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
  void compute(int, int) override;         //Workhorse comuptation
  void settings(int, char **) override;    // Reads args from "pair_style"
  void coeff(int, char **) override;       // Reads args from "pair_coeff" (only * * for mtp)
  void init_style() override;              //Init style
  double init_one(int, int) override;      // Checks that species are inited

 protected:
  void read_file(FILE *);                     //Parsing file using LAMMPS utils
  std::string potential_name = "Untitled";    //An optional name which isn't currently used.
  std::string potential_tag = "";    //An optional tag/description which isn't currently used.

  int species_count;     // Number of species
  double scaling = 1;    // All forces are multiplied by scaling

  // Radial basis
  //1 => "RBChebyshev"
  int radial_basis_type_index;        // Index for MPI Bcast
  RadialMTPBasis *radial_basis;       // Pointer to basis object
  double *radial_basis_coeffs;        // These are the radial basis coeffs (c)
  int radial_func_count;              // Number of radial bases (mu_max)
  int radial_basis_size;              // Number of elements in bases
  int radial_coeff_count;             // Number of total radial coeffs
  int radial_coeff_count_per_pair;    // Number of coeffs for species pair
  double min_cutoff;                  // Min radial cutoff
  double max_cutoff;                  // Max radial cutoff
  double
      max_cutoff_sq;    // Maximum radial cutoff squared (The MTP only supports one cutoff for all species combinations)

  double *linear_coeffs;     // These are the moment tensor basis coeffs (xi)
  double *species_coeffs;    // For the species coefficients (0th rank moment tensor)
  int alpha_moment_count, alpha_index_basic_count, alpha_index_times_count, alpha_scalar_count,
      max_alpha_index_basic;    // Counts of various alpha indicies
  int **alpha_index_basic;      // Indicies how to construct elementary moments from coords and dist
  int **alpha_index_times;      // Indicies to combine existing moments into knew ones
  int *alpha_moment_mapping;    // Selects the basis values from completed moments

  // Other working buffers
  int jac_size = 0;         // Size of the jacobian (jnum dim)
  double *dist_powers;      // Buffer used for powers of dist (eg. d^i)
  double **coord_powers;    // Buffer used for powers of rel. pos. (eg. [dx^i, dy^i, dz^i])
  double *radial_vals;      // Buffer used for radial basis function values for each mu
  double *radial_ders;      // Buffer used for radial basis function derivatives for each mu
  double ***moment_jacobian = nullptr;    // First created during compute using grow
  double *moment_tensor_vals;             //Buffer to hold the moments
  double *nbh_energy_ders_wrt_moments;    // Same as above except for ders

  // Cache whether to calculate forces based on cutoff as calculated in alpha basics
  bool *within_cutoff = nullptr;    // First created during compute using grow
};

}    // namespace LAMMPS_NS

#endif
#endif

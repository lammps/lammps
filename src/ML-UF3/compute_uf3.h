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
ComputeStyle(uf3, ComputeUF3);
// clang-format on
#else

#ifndef LMP_COMPUTE_UF3_H
#define LMP_COMPUTE_UF3_H

#include <string>

#include "compute.h"
#include "pair_uf3.h"

namespace LAMMPS_NS {

/** Minimal PairUF3 wrapper for UF3 compute (2-body descriptor pipeline). */
class UF3Reader : public PairUF3 {
 public:
  explicit UF3Reader(LAMMPS *lmp);
  void load_settings_2body();
  void load_coeffs(int narg, char **argv);
  void ensure_bsplines();
  bool has_3b_pot() const { return pot_3b; }
  double cut_ij(int i, int j) const { return cut[i][j]; }
  int n2b_coeff_count(int i, int j) const { return n2b_coeff_array_size[i][j]; }
  double get_cutsq(int i, int j) const;
  int knot_start_2b(int it, int jt, double r);
  int ncoeff_2b(int i, int j) const;
  double basis_2b_phi(int it, int jt, int lcoeff, double rsq, double rij);
  double basis_2b_dphi_dr(int it, int jt, int lcoeff, double rsq, double rij);
};

class ComputeUF3 : public Compute {
 public:
  ComputeUF3(class LAMMPS *, int, char **);
  ~ComputeUF3() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  double memory_usage() override;

 protected:
  void free_local();

  UF3Reader *uf3_reader;
  double cutmax;
  int lastcol;
  int ncoeff;
  int **col_off_2b;
  class NeighList *list;
  class Compute *c_pe, *c_virial;
  std::string id_virial;
  double **uf3local;
  double **uf3all;

  void build_column_offsets();
};

}    // namespace LAMMPS_NS

#endif
#endif

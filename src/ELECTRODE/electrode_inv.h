/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#ifndef LMP_ELECTRODE_INV_H
#define LMP_ELECTRODE_INV_H

#include "charge_solver.h"
#include "fix.h"
#include "pointers.h"
#include <unordered_map>

namespace LAMMPS_NS {

class ElectrodeInv : public Pointers, public ChargeSolver {
 public:
  // ChargeSolver methods
  ElectrodeInv(class LAMMPS *);
  ~ElectrodeInv() noexcept;    // TODO why do we need noexcept here
  void update_solver(std::vector<tagint>, std::vector<int>) override;
  void set_elyt_pot(double *) override;
  std::vector<double> solve(std::vector<double>) override;
  std::vector<double> compute_potentials() override;
  double get_potential(int) override;
  double get_sb_charges(int) override;
  double get_macro_capacitance(int, int) override;
  double get_macro_elastance(int, int) override;
  void buffer_and_gather(double const *, double *) override;
  double vacuum_capacitance() override;    // for electrode/thermo
  double memory_use() override;

  // setup
  void set_capacitance(int, double **);
  void set_elastance(int, double **);
  void setup_solver(int, std::unordered_map<tagint, int>, std::vector<int>, bool);

 private:
  int groupbit;
  int nmax, nlocalele;
  int ngroups, nele_world;
  bigint elyt_step;
  bool setup, cap_set, vac_cap_computed;
  double evscale, macro_capacitance_sum, vac_cap;
  double *potential_i, *potential_iele, *buf_gathered;
  double **capacitance;
  std::vector<double> qvec;
  std::vector<double> sb_charges;    // group charges w/o potential
  std::vector<double> group_pot;     // group potentials, set during last solve
  std::vector<std::vector<double>> macro_capacitance;
  std::vector<std::vector<double>> macro_elastance;

  int *recvcounts, *displs;    // for MPI-building of iele_gathered

  // pre compute
  std::vector<std::vector<double>> sd_vectors;
  void compute_sd_vectors(), compute_sd_vectors_ffield(std::vector<int>);
  void symmetrize();
  int get_top_group(std::vector<int>);
  void compute_macro_matrices(bool);

  // permanent lists
  std::vector<int> iele_to_group;
  std::unordered_map<tagint, int> tag_to_iele;    // inverse of global taglist:

  // non-permanent lists
  std::vector<tagint> taglist_local;
  std::vector<int> iele_local;     // electrode IDs owned by me
  int *iele_gathered;              // MPIgathered iele_local: all electrode IDs, nproc-ordered
  std::vector<double> buf_iele;    // buffer for electrode properties ordered by iele_local

  // updating methods
  std::vector<double> apply_constraint(std::vector<double>);
};

}    // namespace LAMMPS_NS

#endif


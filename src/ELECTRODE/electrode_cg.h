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

#ifndef LMP_ELECTRODE_CG_H
#define LMP_ELECTRODE_CG_H

#include "charge_solver.h"
#include "electrode_vector.h"
#include "fix.h"
class ElectrodeVector;

namespace LAMMPS_NS {

class ElectrodeCG : public Fix, public ChargeSolver {
 public:
  // ChargeSolver methods
  ElectrodeCG(class LAMMPS *);
  ~ElectrodeCG() noexcept;    // TODO why do we need noexcept here
  void update_solver(std::vector<tagint>, std::vector<int>) override;
  void set_elyt_pot(double *) override;
  std::vector<double> solve(std::vector<double>) override;
  std::vector<double> compute_potentials() override;
  double get_potential(int) override;
  double get_sb_charges(int) override;
  double get_macro_capacitance(int, int) override;
  double get_macro_elastance(int, int) override;
  void buffer_and_gather(double const *, double *) override;
  double memory_use() override;

  // for electrode/thermo (not implemented yet)
  double vacuum_capacitance() override;

  //setup
  void setup_solver(double, ElectrodeVector *, int);

  // fix methods
  int setmask() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  int nele, nele_world;
  virtual void setup_cg(double, int);
  virtual std::vector<double> ele_ele_interaction(const std::vector<double> &);

 private:
  int nmax;
  long nstep, ncall;
  bigint elyt_step;
  bool setup, a_cached_flag;
  double evscale, threshold;
  ElectrodeVector *elec_vec;
  std::vector<double> q_ele;
  int predictor_index, predictor_cols, predictor_count;
  std::vector<std::vector<double>> predictor_weights;
  std::vector<tagint> taglist;
  std::vector<int> iele_to_group;
  double *potential_i;    // potentials, i-indexed (0 for non-electrode atoms)
  std::vector<double> bvec, a_cached;

  void predict_q();
  std::vector<double> pot_to_vector(double *);
  void set_charges(std::vector<double>);

  // math operations with vectors
  std::vector<double> scale_vector(double, std::vector<double>);
  std::vector<double> add(std::vector<double>, std::vector<double>);
  double dot_product(std::vector<double>, std::vector<double>);
  //
  std::vector<double> constraint_projection(std::vector<double>, bool);
};

}    // namespace LAMMPS_NS

#endif


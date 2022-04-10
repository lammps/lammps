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

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

// clang-format off
FixStyle(electrode/conp, FixElectrodeConp);
// clang-format on

#else

#ifndef LMP_FIX_ELECTRODE_CONP_H
#define LMP_FIX_ELECTRODE_CONP_H

#include "electrode_accel_interface.h"
#include "fix.h"

#include <algorithm>
#include <map>

namespace LAMMPS_NS {

class FixElectrodeConp : public Fix {
 public:
  FixElectrodeConp(class LAMMPS *, int, char **);
  ~FixElectrodeConp() override;
  int setmask() override;
  void setup_post_neighbor() override;
  void setup_pre_reverse(int, int) override;
  void pre_force(int) override;
  void pre_reverse(int, int) override;
  double compute_scalar() override;
  double compute_vector(int) override;
  int modify_param(int, char **) override;
  int modify_param(const std::string &);
  void init() override;
  void init_list(int, NeighList *) override;
  void post_constructor() override;    // used by ffield to set up fix efield

 protected:
  virtual void update_psi();
  virtual void pre_update(){};
  virtual void compute_macro_matrices();
  std::vector<double> group_psi;
  std::vector<int> group_bits;
  int num_of_groups;
  bigint ngroup;
  std::vector<int> mpos_to_group;
  std::vector<std::vector<double>> sd_vectors;
  std::vector<double> sb_charges;
  std::vector<int> group_psi_var_ids, group_psi_var_styles;
  std::vector<std::string> group_psi_var_names;
  bool symm;                                           // symmetrize elastance for charge neutrality
  std::vector<std::vector<double>> macro_elastance;    // used by conq
  std::vector<std::vector<double>> macro_capacitance;    // used by thermo
  double thermo_temp, thermo_time;                       // used by electrode/thermo only
  bool ffield;                                           // possibly tweak electrode/conq's version
  std::string fixname;    // used by electrode/ffield to set up internal efield
  bool intelflag;
  ElectrodeAccelInterface *accel_interface;    // used by /intel

 private:
  FILE *f_inv, *f_mat, *f_vec;    // files for capacitance, eleastance and vector
  std::string input_file_inv, input_file_mat;
  class ElectrodeMatrix *array_compute;
  class ElectrodeVector *ele_vector;
  std::vector<int> groups;
  double **capacitance;
  std::vector<tagint> taglist, taglist_bygroup, group_idx;
  std::map<tagint, int> tag_to_iele;
  bool read_inv, read_mat;
  double eta;
  double update_time, mult_time;
  void create_taglist();
  void invert();
  void symmetrize();
  double gausscorr(int, bool);
  void update_charges();
  double potential_energy(int, const std::vector<int> &);
  double self_energy(int);
  std::vector<int> local_to_matrix();
  void write_to_file(FILE *, const std::vector<tagint> &, const std::vector<std::vector<double>> &);
  void read_from_file(std::string input_file, double **, const std::string &);
  void compute_sd_vectors();
  void compute_sd_vectors_ffield();
  std::vector<int> setvars_types, setvars_groups, setvars_vars;
  void update_setvars(int);
  int groupnum_from_name(char *);
  double evscale;
  class Pair *pair;
  class NeighList *mat_neighlist, *vec_neighlist;
  std::vector<int> etypes;
  int mat_request, vec_request;
  void request_etypes_neighlists();
  bool etypes_neighlists;
  int get_top_group();    // used by ffield
  int top_group;          // used by ffield
  bool tfflag;
  bool timer_flag;
  std::map<int, double> tf_types;
};

}    // namespace LAMMPS_NS

#endif
#endif

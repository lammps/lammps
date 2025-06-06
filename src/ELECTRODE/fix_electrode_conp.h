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

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (GU), Kamila Savvidi (TUHH), Robert Meissner (Hereon, TUHH)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

// clang-format off
FixStyle(electrode/conp, FixElectrodeConp);
// clang-format on

#else

#ifndef LMP_FIX_ELECTRODE_CONP_H
#define LMP_FIX_ELECTRODE_CONP_H

#include "fix.h"

#include <map>

namespace LAMMPS_NS {
// forward decls

class ChargeSolver;
class ElectrodeTaglist;
class ElectrodeVector;
class NeighList;
class Pair;

class FixElectrodeConp : public Fix {

 public:
  FixElectrodeConp(class LAMMPS *, int, char **);
  ~FixElectrodeConp() override;
  int setmask() override;
  void setup_pre_exchange() override;
  void setup_post_neighbor() override;
  void setup_pre_reverse(int, int) override;
  void pre_force(int) override;
  void pre_reverse(int, int) override;
  double compute_scalar() override;
  double compute_vector(int) override;
  double compute_array(int, int) override;
  int modify_param(int, char **) override;
  int modify_param(const std::string &);
  virtual void init() override;
  void init_list(int, NeighList *) override;
  void post_constructor() override;    // used by ffield to set up fix efield
  double memory_usage() override;

  // atomvec-based tracking of electrode atoms
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  enum class Algo { MATRIX_INV, MATRIX_CG, CG };
  enum class VarStyle { CONST, EQUAL, UNSET };
  ChargeSolver *charge_solver;
  virtual void update_psi_set_constraint();
  std::vector<double> group_psi;
  std::vector<double> group_psi_const;    // needed to undo qtotal psi updates
  std::vector<int> group_bits;
  std::vector<int> groups;
  int num_of_groups;
  bigint ngroup;
  double evscale;
  double *potential_i;    // potentials, i-indexed (0 for non-electrode atoms)
  std::vector<int> group_psi_var_ids;
  std::vector<VarStyle> group_psi_var_styles;
  std::vector<std::string> group_psi_var_names;
  Algo algo;
  double thermo_temp, thermo_time;    // used by electrode/thermo only
  int thermo_init;                    // initializer for rng in electrode/thermo
  bool ffield;                        // possibly tweak electrode/conq's version
  std::string fixname;                // used by electrode/ffield to set up internal efield
  bool intelflag;
  inline virtual void intel_pack_buffers() {}
  double qtotal;
  std::string qtotal_var_name;
  int qtotal_var_id;
  VarStyle qtotal_var_style;

 private:
  std::string output_file_inv, output_file_mat, output_file_vec;
  std::string input_file_inv, input_file_mat;
  ElectrodeVector *elyt_vector, *elec_vector;
  double **matrix;
  bool read_inv, read_mat, write_inv, write_mat, write_vec;
  bool matrix_algo, need_array_compute, need_elec_vector;
  double eta, cg_threshold;
  double update_time, mult_time;
  double gausscorr(int, int, bool);
  void update_charges();
  double potential_energy();
  double self_energy(int);
  void v_tally(int, int, int, int, double, double, double, double);
  Pair *pair;
  NeighList *mat_neighlist, *vec_neighlist;
  std::vector<int> etypes;
  void request_etypes_neighlists();
  bool etypes_neighlists;
  int get_top_group();    // used by ffield
  int top_group;          // used by ffield
  bool tfflag;
  void add_electronegativity(double *);
  bool enflag, hardnessflag;                  // qeq parameters set
  int eta_index, hardness_index, en_index;    // index of atom property for eta
  bool etapropflag;                           // eta specified as atom property
  bool pairflag;                              // whether a pair style is specified
  std::string pair_str;
  bool timer_flag;
  std::map<int, double> tf_types;
  void set_charges(std::vector<double>);

  // fix-specific electrode ID storage system:
  bool taglist_constructed;
  ElectrodeTaglist *electrode_taglist;
  int nlocalele_outdated;    // trigger rebuilding of following structures:
  std::vector<tagint> taglist_local;
  std::vector<int> iele_to_group_local;

  int nlocalele;              // current no. of local electrode atoms
  void gather_list_iele();    // build iele_gathered

  int nmax;
};

}    // namespace LAMMPS_NS

#endif
#endif

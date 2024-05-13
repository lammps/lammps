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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
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
#include <unordered_map>

namespace LAMMPS_NS {
// forward decls

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

  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  enum class Algo { MATRIX_INV, MATRIX_CG, CG };
  enum class VarStyle { CONST, EQUAL, UNSET };
  virtual void update_psi();
  virtual void pre_update(){};
  virtual void recompute_potential(std::vector<double>, std::vector<double>){};
  virtual std::vector<double> constraint_projection(std::vector<double>);
  virtual std::vector<double> constraint_correction(std::vector<double>);
  virtual void compute_macro_matrices();
  std::vector<double> ele_ele_interaction(const std::vector<double> &);
  std::vector<double> group_psi;
  std::vector<double> group_psi_const;    // needed to undo qtotal psi updates
  std::vector<int> group_bits;
  std::vector<int> groups;
  int num_of_groups;
  bigint ngroup;
  double evscale;
  std::vector<std::vector<double>> sd_vectors;
  std::vector<double> sb_charges;
  std::vector<int> group_psi_var_ids;
  std::vector<VarStyle> group_psi_var_styles;
  std::vector<std::string> group_psi_var_names;
  std::vector<int> iele_to_group_local;
  bool symm;    // symmetrize elastance for charge neutrality
  Algo algo;
  std::vector<std::vector<double>> macro_elastance;      // used by conq
  std::vector<std::vector<double>> macro_capacitance;    // used by thermo
  double thermo_temp, thermo_time;                       // used by electrode/thermo only
  int thermo_init;                                       // initializer for rng in electrode/thermo
  bool ffield;                                           // possibly tweak electrode/conq's version
  std::string fixname;    // used by electrode/ffield to set up internal efield
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
  double **capacitance, **elastance;
  bool read_inv, read_mat, write_inv, write_mat, write_vec;
  bool matrix_algo, need_array_compute, need_elec_vector;
  double eta, cg_threshold;
  double update_time, mult_time;
  long n_cg_step, n_call;
  void create_taglist();
  void invert();
  void symmetrize();
  double gausscorr(int, bool);
  void update_charges();
  double potential_energy();
  double self_energy(int);
  void write_to_file(FILE *, const std::vector<tagint> &, const std::vector<std::vector<double>> &);
  void read_from_file(const std::string &input_file, double **, const std::string &);
  void compute_sd_vectors();
  void compute_sd_vectors_ffield();
  int groupnum_from_name(char *);
  Pair *pair;
  NeighList *mat_neighlist, *vec_neighlist;
  std::vector<int> etypes;
  int mat_request, vec_request;
  void request_etypes_neighlists();
  bool etypes_neighlists;
  int get_top_group();    // used by ffield
  int top_group;          // used by ffield
  bool tfflag;
  int eta_index;    // index of atom property for eta
  bool etaflag;     // eta specified as atom property
  bool timer_flag;
  std::map<int, double> tf_types;
  // cg
  std::vector<double> scale_vector(double, std::vector<double>);
  std::vector<double> add_nlocalele(std::vector<double>, std::vector<double>);
  double dot_nlocalele(std::vector<double>, std::vector<double>);
  std::vector<double> times_elastance(std::vector<double>);
  std::vector<double> gather_ngroup(std::vector<double>);
  std::vector<double> gather_elevec_local(ElectrodeVector *);
  void set_charges(std::vector<double>);
  // qtotal
  double macro_capacitance_sum;
  void update_psi_qtotal();

  // fix-specific electrode ID storage system:

  std::vector<tagint> taglist;    // global list: all tags in combined electrode group
  std::vector<tagint> taglist_local;
  std::vector<tagint> taglist_bygroup;            // taglist sorted by group
  std::vector<tagint> group_idx;                  // permutation taglist<->taglist_bygroup
  std::unordered_map<tagint, int> tag_to_iele;    // inverse of taglist:
  std::vector<int> iele_to_group;
  // tag_to_iele[taglist[iele]] = iele
  int nlocalele;    // current no. of local electrode atoms

  int nlocalele_outdated;          // trigger rebuilding of following structures:
  std::vector<int> list_iele;      // electrode IDs owned by me
  int *recvcounts, *displs;        // for MPI-building of iele_gathered
  int *iele_gathered;              // MPIgathered list_iele: all electrode IDs, nproc-ordered
  std::vector<double> buf_iele;    // buffer for electrode properties ordered by list_iele
  double *buf_gathered;            // buffer for MPIgathered buf_iele (NOT YET iele-ordered)
  double *potential_i;             // potentials, i-indexed (0 for non-electrode atoms)
  double *potential_iele;          // potentials ordered by iele

  void gather_list_iele();                       // build iele_gathered
  void gather_elevec(double *);                  // gather buf_iele and rearrange into iele-order
  void buffer_and_gather(double *, double *);    // buffer into buf_iele then gather and rearrange

  int nmax;
};

}    // namespace LAMMPS_NS

#endif
#endif

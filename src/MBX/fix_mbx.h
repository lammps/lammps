/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(MBX, FixMBX)

#else

#ifndef LMP_FIX_MBX_H
#define LMP_FIX_MBX_H

#include "fix.h"



namespace LAMMPS_NS {

class FixMBX : public Fix {
  friend class PairMBX;

 public:
  FixMBX(class LAMMPS *, int, char **);
  ~FixMBX();
  int setmask();
  virtual void post_constructor();
  virtual void init();
  virtual void init_storage();
  void setup(int);
  void min_setup(int);
  virtual void setup_post_neighbor();
  virtual void post_neighbor();
  void min_post_neighbor();
  void setup_pre_force(int);
  virtual void pre_force(int);

  void min_setup_pre_force(int);
  void min_pre_force(int);

  void setup_pre_exchange();
  void pre_exchange();

  void post_force(int);
  void min_post_force(int);

 protected:
  struct MBXImpl *mbx_impl;
  class PairMBX *pair_mbx;    // pointer to MBX pair_style

  static std::string cite_pair_mbx;
  int me, nprocs;
  bigint ngroup;

  bool mbx_aspc_enabled;
  bool print_dipoles;

  bool first_step;
  bool has_gcmc;

  int use_json;
  std::string json_settings;

  int print_verbose;

  int num_mol_types;    // # of unique molecule types
  int num_molecules;    // total # of molecules
  int *num_atoms_per_mol;                // array of # of atoms per molecule for each type
  int *lower_atom_type_index_in_mol;     // array with the lowest atom type index in the monomer
  int *higher_atom_type_index_in_mol;    // array with the highest atom type index in the monomer
  int **order_in_mol;                    // array with the atom order for each monomer
  char **mol_names;    // array of molecule names

  int *mol_type;      // per-atom array of molecule type
  int *mol_anchor;    // per-atom array 1/0 if anchor atom of a molecule
  int *mol_local;    // per-molecule array 1/0 if molecule has at least one local particle

  int mbx_num_atoms, mbx_num_ext;
  int mbx_num_atoms_local, mbx_num_ext_local;

  int *mbxt_count;
  double *mbxt_time;
  double *mbxt_time_start;
  double mbxt_initial_time;
  void mbxt_start(int);
  void mbxt_stop(int);
  void mbxt_write_summary();
  void mbxt_print_time(const char *, int, double *);
  enum MBXT_LABELS : int {
    INIT = 0,
    UPDATE_XYZ,
    INIT_LOCAL,
    UPDATE_XYZ_LOCAL,
    E1B,
    E2B_GHOST,
    E3B_GHOST,
    E4B_GHOST,
    DISP,
    DISP_PME,
    BUCK,
    ELE,
    ACCUMULATE_F,
    ACCUMULATE_F_LOCAL,

    ELE_PERMDIP_REAL,
    ELE_PERMDIP_PME,

    ELE_DIPFIELD_REAL,
    ELE_DIPFIELD_PME,

    ELE_GRAD_REAL,
    ELE_GRAD_PME,
    ELE_GRAD_FIN,

    ELE_COMM_REVFOR,
    ELE_COMM_REVSET,
    ELE_COMM_REV,
    ELE_COMM_FORSET,
    ELE_COMM_FOR,

    ELE_PME_SETUP,
    ELE_PME_C,
    ELE_PME_D,
    ELE_PME_E,

    DISP_PME_SETUP,
    DISP_PME_E,

    NUM_TIMERS
  };

  bool validateMBXFixParameters(int narg, char **arg);

  int aspc_order;
  int aspc_num_hist;
  int aspc_max_num_hist;
  int aspc_per_atom_size;
  double **aspc_dip_hist;

  int aspc_step;
  int aspc_step_reset;

  double **mbx_dip;

  void mbx_init();
  void mbx_init_local();

  void mbx_update_xyz();
  void mbx_update_xyz_local();

  void mbx_fill_system_information_from_atom();

  void mbx_init_dipole_history_local();

  void mbx_get_dipoles_local();

  int get_num_atoms_per_monomer(char *, bool &);
  int get_include_monomer(char *, int, bool &, bool &);
  void add_monomer_atom_types(char *, std::vector<std::string> &);
  std::pair<int, int> parse_dp1_range(const std::string &);

  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  virtual void grow_arrays(int);
  virtual void copy_arrays(int, int, int);
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(int, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif

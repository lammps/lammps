/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Heavily modified and adapted for LAMMPS by the LAMMPS developers.
------------------------------------------------------------------------- */

#ifndef LMP_REAXFF_TYPES_H
#define LMP_REAXFF_TYPES_H

#include "lmptype.h"

#include "reaxff_defs.h"      // IWYU pragma: export
#include "reaxff_inline.h"    // IWYU pragma: export

// forward declarations
namespace LAMMPS_NS {
class Error;
class LAMMPS;
class Memory;
class Pair;
}    // namespace LAMMPS_NS

namespace ReaxFF {
/********************** TYPE DEFINITIONS ********************/
typedef int ivec[3];
typedef double rvec[3];
typedef double rvec2[2];

// import LAMMPS' definition of tagint and bigint
typedef LAMMPS_NS::tagint rc_tagint;
typedef LAMMPS_NS::bigint rc_bigint;

struct global_parameters {
  int n_global;
  int vdw_type;
  double *l;
};

struct single_body_parameters {
  char name[4];    // two character atom name
  double r_s;
  double valency;    // Valency of the atom
  double mass;       // Mass of atom
  double r_vdw;
  double epsilon;
  double gamma;
  double r_pi;
  double valency_e;
  double nlp_opt;

  /* Line two in field file */
  double alpha;
  double gamma_w;
  double valency_boc;
  double p_ovun5;
  double chi;
  double eta;
  int p_hbond;    // 1 for H, 2 for hbonding atoms (O,S,P,N), 0 for others

  /* Line three in field file */
  double r_pi_pi;
  double p_lp2;
  double b_o_131;
  double b_o_132;
  double b_o_133;
  double bcut_acks2;    // ACKS2 bond cutoff

  /* Line four in the field file */
  double p_ovun2;
  double p_val3;
  double valency_val;
  double p_val5;
  double rcore2;
  double ecore2;
  double acore2;

  /* Line five in the ffield file, only for lgvdw yes */
  double lgcij;
  double lgre;
};

/* Two Body Parameters */
struct two_body_parameters {
  /* Bond Order parameters */
  double p_bo1, p_bo2, p_bo3, p_bo4, p_bo5, p_bo6;
  double r_s, r_p, r_pp;    // r_o distances in BO formula
  double p_boc3, p_boc4, p_boc5;

  /* Bond Energy parameters */
  double p_be1, p_be2;
  double De_s, De_p, De_pp;

  /* Over/Under coordination parameters */
  double p_ovun1;

  /* Van der Waal interaction parameters */
  double D;
  double alpha;
  double r_vdW;
  double gamma_w;
  double rcore, ecore, acore;
  double lgcij, lgre;

  /* electrostatic parameters */
  double gamma;    // note: this parameter is gamma^-3 and not gamma.

  double v13cor, ovc;
};

struct dbond_coefficients {
  double C1dbo, C2dbo, C3dbo;
  double C1dbopi, C2dbopi, C3dbopi, C4dbopi;
  double C1dbopi2, C2dbopi2, C3dbopi2, C4dbopi2;
  double C1dDelta, C2dDelta, C3dDelta;
};

/* 3-body parameters */
struct three_body_parameters {
  /* valence angle */
  double theta_00;
  double p_val1, p_val2, p_val4, p_val7;

  /* penalty */
  double p_pen1;

  /* 3-body conjugation */
  double p_coa1;
};

struct three_body_header {
  int cnt;
  three_body_parameters prm[REAX_MAX_3BODY_PARAM];
};

/* hydrogen-bond parameters */
struct hbond_parameters {
  double r0_hb, p_hb1, p_hb2, p_hb3;
};

/* 4-body parameters */
struct four_body_parameters {
  double V1, V2, V3;

  /* torsion angle */
  double p_tor1;

  /* 4-body conjugation */
  double p_cot1;
};

struct four_body_header {
  int cnt;
  four_body_parameters prm[REAX_MAX_4BODY_PARAM];
};

struct reax_interaction {
  int num_atom_types;
  global_parameters gp;
  single_body_parameters *sbp;
  two_body_parameters **tbp;
  three_body_header ***thbp;
  hbond_parameters ***hbp;
  four_body_header ****fbp;
};

struct reax_atom {
  rc_tagint orig_id;
  int type;
  char name[8];

  rvec x;      // position
  rvec v;      // velocity
  rvec f;      // force
  double q;    // charge

  int Hindex;
  int num_bonds;
  int num_hbonds;
};

struct LR_lookup_table;    // forward declaration
struct reax_system {
  reax_interaction reax_param;

  int n, N, numH;
  int local_cap, total_cap, Hcap;
  int wsize, my_rank, num_nbrs;
  reax_atom *my_atoms;

  LAMMPS_NS::Error *error_ptr;
  LAMMPS_NS::Pair *pair_ptr;
  LAMMPS_NS::Memory *mem_ptr;

  int my_bonds;
  int mincap, minhbonds;
  double safezone, saferzone;

  LR_lookup_table **LR;

  int omp_active;
  int acks2_flag;
};

/* system control parameters */
struct control_params {
  int nthreads;

  double bond_cut;
  double nonb_cut, nonb_low;
  double hbond_cut;

  double bg_cut;
  double bo_cut;
  double thb_cut;
  double thb_cutsq;

  int tabulate;

  int lgflag;
  int enobondsflag;
  LAMMPS_NS::Error *error_ptr;
  LAMMPS_NS::LAMMPS *lmp_ptr;
  int me;
};

struct energy_data {
  double e_bond;    // Total bond energy
  double e_ov;      // Total over coordination
  double e_un;      // Total under coordination energy
  double e_lp;      // Total under coordination energy
  double e_ang;     // Total valance angle energy
  double e_pen;     // Total penalty energy
  double e_coa;     // Total three body conjgation energy
  double e_hb;      // Total Hydrogen bond energy
  double e_tor;     // Total torsional energy
  double e_con;     // Total four body conjugation energy
  double e_vdW;     // Total van der Waals energy
  double e_ele;     // Total electrostatics energy
  double e_pol;     // Polarization energy
};

struct simulation_data {
  rc_bigint step;
  energy_data my_en;    // per MPI rank energies
};

struct three_body_interaction_data {
  int thb;
  int pthb;    // pointer to the third body on the central atom's nbrlist
  double theta, cos_theta;
  rvec dcos_di, dcos_dj, dcos_dk;
};

struct far_neighbor_data {
  int nbr;
  ivec rel_box;
  double d;
  rvec dvec;
};

struct hbond_data {
  int nbr;
  int scl;
  far_neighbor_data *ptr;
};

struct bond_order_data {
  double BO, BO_s, BO_pi, BO_pi2;
  double Cdbo, Cdbopi, Cdbopi2;
  double C1dbo, C2dbo, C3dbo;
  double C1dbopi, C2dbopi, C3dbopi, C4dbopi;
  double C1dbopi2, C2dbopi2, C3dbopi2, C4dbopi2;
  rvec dBOp, dln_BOp_s, dln_BOp_pi, dln_BOp_pi2;
  double *CdboReduction;
};

struct bond_data {
  int nbr;
  int sym_index;
  int dbond_index;
  ivec rel_box;
  double d;
  rvec dvec;
  bond_order_data bo_data;
};

struct sparse_matrix_entry {
  int j;
  double val;
};

struct sparse_matrix {
  int cap, n, m;
  int *start, *end;
  sparse_matrix_entry *entries;
};

struct reallocate_data {
  int num_far;
  int H, Htop;
  int hbonds, num_hbonds;
  int bonds, num_bonds;
  int num_3body;
};

struct storage {
  int allocated;

  /* bond order related storage */
  double *total_bond_order;
  double *Deltap, *Deltap_boc;
  double *Delta, *Delta_lp, *Delta_lp_temp, *Delta_e, *Delta_boc, *Delta_val;
  double *dDelta_lp, *dDelta_lp_temp;
  double *nlp, *nlp_temp, *Clp, *vlpex;
  rvec *dDeltap_self;
  int *bond_mark;

  /* Taper */
  double Tap[8];    //Tap7, Tap6, Tap5, Tap4, Tap3, Tap2, Tap1, Tap0;

  /* force calculations */
  double *CdDelta;    // coefficient of dDelta
  rvec *f;

  /* omp */
  rvec *forceReduction;
  double *CdDeltaReduction;
  int *valence_angle_atom_myoffset;

  /* acks2 */
  double *s;

  reallocate_data realloc;
};

union list_type {
  three_body_interaction_data *three_body_list;
  bond_data *bond_list;
  far_neighbor_data *far_nbr_list;
  hbond_data *hbond_list;
};

struct reax_list {
  int allocated;

  int n;
  int num_intrs;

  int *index;
  int *end_index;

  int type;
  list_type select;
  class LAMMPS_NS::Error *error_ptr;
};

struct LR_lookup_table {
  double xmin, xmax;
  int n;
  double dx, inv_dx;
  double a;
  double m;
  double c;

  LR_data *y;
  cubic_spline_coef *H;
  cubic_spline_coef *vdW, *CEvd;
  cubic_spline_coef *ele, *CEclmb;
};
}    // namespace ReaxFF

#endif

/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#ifndef __REAX_TYPES_H_
#define __REAX_TYPES_H_

#include <mpi.h>
#include "lmptype.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/time.h>
#include "accelerator_kokkos.h"


namespace LAMMPS_NS { class Error;}

#if defined LMP_USER_OMP
#define OMP_TIMING 0

#ifdef OMP_TIMING
// pkcoff timing fields
enum {
        COMPUTEINDEX=0,
        COMPUTEWLINDEX,
        COMPUTEBFINDEX,
        COMPUTEQEQINDEX,
        COMPUTENBFINDEX,
        COMPUTEIFINDEX,
        COMPUTETFINDEX,
        COMPUTEBOINDEX,
        COMPUTEBONDSINDEX,
        COMPUTEATOMENERGYINDEX,
        COMPUTEVALENCEANGLESBOINDEX,
        COMPUTETORSIONANGLESBOINDEX,
        COMPUTEHBONDSINDEX,
        COMPUTECG1INDEX,
        COMPUTECG2INDEX,
        COMPUTECGCOMPUTEINDEX,
        COMPUTECALCQINDEX,
        COMPUTEINITMVINDEX,
        COMPUTEMVCOMPINDEX,
        LASTTIMINGINDEX
};

extern double ompTimingData[LASTTIMINGINDEX];
extern int ompTimingCount[LASTTIMINGINDEX];
extern int ompTimingCGCount[LASTTIMINGINDEX];
#endif
#endif

/************* SOME DEFS - crucial for reax_types.h *********/

#define LAMMPS_REAX

//#define DEBUG
//#define DEBUG_FOCUS
//#define TEST_ENERGY
//#define TEST_FORCES
//#define CG_PERFORMANCE
//#define LOG_PERFORMANCE
//#define STANDARD_BOUNDARIES
//#define OLD_BOUNDARIES
//#define MIDPOINT_BOUNDARIES

#define REAX_MAX_STR            1024
#define REAX_MAX_NBRS           6
#define REAX_MAX_3BODY_PARAM    5
#define REAX_MAX_4BODY_PARAM    5
#define REAX_MAX_ATOM_TYPES     25
#define REAX_MAX_MOLECULE_SIZE  20
#define MAX_BOND                    20  // same as reaxc_defs.h

/********************** TYPE DEFINITIONS ********************/
typedef int ivec[3];
typedef double rvec[3];
typedef double rtensor[3][3];
typedef double rvec2[2];
typedef double rvec4[4];

// import LAMMPS' definition of tagint and bigint
typedef LAMMPS_NS::tagint rc_tagint;
typedef LAMMPS_NS::bigint rc_bigint;

typedef struct
{
  int  cnt;
  int *index;
  void *out_atoms;
} mpi_out_data;

typedef struct
{
  MPI_Comm     world;
  MPI_Comm     comm_mesh3D;

  MPI_Datatype sys_info;
  MPI_Datatype mpi_atom_type;
  MPI_Datatype boundary_atom_type;
  MPI_Datatype mpi_rvec, mpi_rvec2;
  MPI_Datatype restart_atom_type;

  MPI_Datatype header_line;
  MPI_Datatype header_view;
  MPI_Datatype init_desc_line;
  MPI_Datatype init_desc_view;
  MPI_Datatype atom_line;
  MPI_Datatype atom_view;
  MPI_Datatype bond_line;
  MPI_Datatype bond_view;
  MPI_Datatype angle_line;
  MPI_Datatype angle_view;

  mpi_out_data out_buffers[REAX_MAX_NBRS];
  void *in1_buffer;
  void *in2_buffer;
} mpi_datatypes;

typedef struct
{
  int n_global;
  double* l;
  int vdw_type;
} global_parameters;

typedef struct
{
  /* Line one in field file */
  char name[15]; // Two character atom name

  double r_s;
  double valency;  // Valency of the atom
  double mass;     // Mass of atom
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
  int  p_hbond; // 1 for H, 2 for hbonding atoms (O,S,P,N), 0 for others

  /* Line three in field file */
  double r_pi_pi;
  double p_lp2;
  double b_o_131;
  double b_o_132;
  double b_o_133;

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

} single_body_parameters;

/* Two Body Parameters */
typedef struct {
  /* Bond Order parameters */
  double p_bo1,p_bo2,p_bo3,p_bo4,p_bo5,p_bo6;
  double r_s, r_p, r_pp;  // r_o distances in BO formula
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
  double gamma; // note: this parameter is gamma^-3 and not gamma.

  double v13cor, ovc;
} two_body_parameters;

/* 3-body parameters */
typedef struct {
  /* valence angle */
  double theta_00;
  double p_val1, p_val2, p_val4, p_val7;

  /* penalty */
  double p_pen1;

  /* 3-body conjugation */
  double p_coa1;
} three_body_parameters;


typedef struct{
  int cnt;
  three_body_parameters prm[REAX_MAX_3BODY_PARAM];
} three_body_header;

/* hydrogen-bond parameters */
typedef struct{
  double r0_hb, p_hb1, p_hb2, p_hb3;
} hbond_parameters;

/* 4-body parameters */
typedef struct {
  double V1, V2, V3;

  /* torsion angle */
  double p_tor1;

  /* 4-body conjugation */
  double p_cot1;
} four_body_parameters;

typedef struct
{
  int cnt;
  four_body_parameters prm[REAX_MAX_4BODY_PARAM];
} four_body_header;

typedef struct
{
  int num_atom_types;
  global_parameters gp;
  single_body_parameters *sbp;
  two_body_parameters **tbp;
  three_body_header ***thbp;
  hbond_parameters ***hbp;
  four_body_header ****fbp;
} reax_interaction;

struct _reax_atom
{
  rc_tagint  orig_id;
  int  imprt_id;
  int  type;
  char name[8];

  rvec x; // position
  rvec v; // velocity
  rvec f; // force
  rvec f_old;

  double q; // charge
  rvec4 s; // they take part in
  rvec4 t; // computing q

  int Hindex;
  int num_bonds;
  int num_hbonds;
  int renumber;

  int numbonds;                  // true number of bonds around atoms
  int nbr_id[MAX_BOND];          // ids of neighbors around atoms
  double nbr_bo[MAX_BOND];          // BO values of bond between i and nbr
  double sum_bo, no_lp;           // sum of BO values and no. of lone pairs
};
typedef _reax_atom reax_atom;

typedef struct
{
  double V;
  rvec min, max, box_norms;

  rtensor box, box_inv;
  rtensor trans, trans_inv;
  rtensor g;
} simulation_box;

struct grid_cell
{
  double cutoff;
  rvec min, max;
  ivec rel_box;

  int  mark;
  int  type;
  int  str;
  int  end;
  int  top;
  int* atoms;
  struct grid_cell** nbrs;
  ivec* nbrs_x;
  rvec* nbrs_cp;
};

typedef struct grid_cell grid_cell;


typedef struct
{
  int  total, max_atoms, max_nbrs;
  ivec ncells;
  rvec cell_len;
  rvec inv_len;

  ivec bond_span;
  ivec nonb_span;
  ivec vlist_span;

  ivec native_cells;
  ivec native_str;
  ivec native_end;

  double ghost_cut;
  ivec ghost_span;
  ivec ghost_nonb_span;
  ivec ghost_hbond_span;
  ivec ghost_bond_span;

  grid_cell*** cells;
  ivec *order;
} grid;


typedef struct
{
  int  rank;
  int  est_send, est_recv;
  int  atoms_str, atoms_cnt;
  ivec rltv, prdc;
  rvec bndry_min, bndry_max;

  int  send_type;
  int  recv_type;
  ivec str_send;
  ivec end_send;
  ivec str_recv;
  ivec end_recv;
} neighbor_proc;



typedef struct
{
  int N;
  int exc_gcells;
  int exc_atoms;
} bound_estimate;



typedef struct
{
  double ghost_nonb;
  double ghost_hbond;
  double ghost_bond;
  double ghost_cutoff;
} boundary_cutoff;

struct _reax_system
{
  reax_interaction reax_param;

  rc_bigint        bigN;
  int              n, N, numH;
  int              local_cap, total_cap, gcell_cap, Hcap;
  int              est_recv, est_trans, max_recved;
  int              wsize, my_rank, num_nbrs;
  ivec             my_coords;
  neighbor_proc    my_nbrs[REAX_MAX_NBRS];
  int             *global_offset;
  simulation_box   big_box, my_box, my_ext_box;
  grid             my_grid;
  boundary_cutoff  bndry_cuts;
  reax_atom       *my_atoms;

  class LAMMPS_NS::Error *error_ptr;
  class LAMMPS_NS::Pair *pair_ptr;
  int my_bonds;
  int mincap;
  double safezone, saferzone;

  int omp_active;
};
typedef _reax_system reax_system;



/* system control parameters */
typedef struct
{
  char sim_name[REAX_MAX_STR];
  int  nprocs;
  int  nthreads;
  ivec procs_by_dim;
  /* ensemble values:
     0 : NVE
     1 : bNVT (Berendsen)
     2 : nhNVT (Nose-Hoover)
     3 : sNPT (Parrinello-Rehman-Nose-Hoover) semiisotropic
     4 : iNPT (Parrinello-Rehman-Nose-Hoover) isotropic
     5 : NPT  (Parrinello-Rehman-Nose-Hoover) Anisotropic*/
  int  ensemble;
  int  nsteps;
  double dt;
  int  geo_format;
  int  restart;

  int  restrict_bonds;
  int  remove_CoM_vel;
  int  random_vel;
  int  reposition_atoms;

  int  reneighbor;
  double vlist_cut;
  double bond_cut;
  double nonb_cut, nonb_low;
  double hbond_cut;
  double user_ghost_cut;

  double bg_cut;
  double bo_cut;
  double thb_cut;
  double thb_cutsq;

  int tabulate;

  int qeq_freq;
  double q_err;
  int refactor;
  double droptol;

  double T_init, T_final, T;
  double Tau_T;
  int  T_mode;
  double T_rate, T_freq;

  int  virial;
  rvec P, Tau_P, Tau_PT;
  int  press_mode;
  double compressibility;

  int  molecular_analysis;
  int  num_ignored;
  int  ignore[REAX_MAX_ATOM_TYPES];

  int  dipole_anal;
  int  freq_dipole_anal;
  int  diffusion_coef;
  int  freq_diffusion_coef;
  int  restrict_type;

  int lgflag;
  int enobondsflag;
  class LAMMPS_NS::Error *error_ptr;
  int me;

} control_params;


typedef struct
{
  double T;
  double xi;
  double v_xi;
  double v_xi_old;
  double G_xi;

} thermostat;


typedef struct
{
  double P;
  double eps;
  double v_eps;
  double v_eps_old;
  double a_eps;

} isotropic_barostat;


typedef struct
{
  rtensor P;
  double P_scalar;

  double eps;
  double v_eps;
  double v_eps_old;
  double a_eps;

  rtensor h0;
  rtensor v_g0;
  rtensor v_g0_old;
  rtensor a_g0;

} flexible_barostat;


typedef struct
{
  double start;
  double end;
  double elapsed;

  double total;
  double comm;
  double nbrs;
  double init_forces;
  double bonded;
  double nonb;
  double qEq;
  int  s_matvecs;
  int  t_matvecs;
} reax_timing;


typedef struct
{
  double e_tot;
  double e_kin;                      // Total kinetic energy
  double e_pot;

  double e_bond;                     // Total bond energy
  double e_ov;                       // Total over coordination
  double e_un;                       // Total under coordination energy
  double e_lp;                       // Total under coordination energy
  double e_ang;                      // Total valance angle energy
  double e_pen;                      // Total penalty energy
  double e_coa;                      // Total three body conjgation energy
  double e_hb;                       // Total Hydrogen bond energy
  double e_tor;                      // Total torsional energy
  double e_con;                      // Total four body conjugation energy
  double e_vdW;                      // Total van der Waals energy
  double e_ele;                      // Total electrostatics energy
  double e_pol;                      // Polarization energy
} energy_data;

typedef struct
{
  int  step;
  int  prev_steps;
  double time;

  double M;                           // Total Mass
  double inv_M;                      // 1 / Total Mass

  rvec xcm;                        // Center of mass
  rvec vcm;                        // Center of mass velocity
  rvec fcm;                        // Center of mass force
  rvec amcm;                       // Angular momentum of CoM
  rvec avcm;                       // Angular velocity of CoM
  double etran_cm;                   // Translational kinetic energy of CoM
  double erot_cm;                    // Rotational kinetic energy of CoM

  rtensor kinetic;                 // Kinetic energy tensor
  rtensor virial;                  // Hydrodynamic virial

  energy_data my_en;
  energy_data sys_en;

  double               N_f;          //Number of degrees of freedom
  rvec               t_scale;
  rtensor            p_scale;
  thermostat         therm;        // Used in Nose_Hoover method
  isotropic_barostat iso_bar;
  flexible_barostat  flex_bar;
  double               inv_W;

  double kin_press;
  rvec int_press;
  rvec my_ext_press;
  rvec ext_press;
  rvec tot_press;

  reax_timing timing;
} simulation_data;


typedef struct{
  int thb;
  int pthb; // pointer to the third body on the central atom's nbrlist
  double theta, cos_theta;
  rvec dcos_di, dcos_dj, dcos_dk;
} three_body_interaction_data;


typedef struct {
  int nbr;
  ivec rel_box;
  double d;
  rvec dvec;
} far_neighbor_data;


typedef struct {
  int nbr;
  int scl;
  far_neighbor_data *ptr;
} hbond_data;


typedef struct{
  int wrt;
  rvec dVal;
} dDelta_data;


typedef struct{
  int wrt;
  rvec dBO, dBOpi, dBOpi2;
} dbond_data;

typedef struct{
  double BO, BO_s, BO_pi, BO_pi2;
  double Cdbo, Cdbopi, Cdbopi2;
  double C1dbo, C2dbo, C3dbo;
  double C1dbopi, C2dbopi, C3dbopi, C4dbopi;
  double C1dbopi2, C2dbopi2, C3dbopi2, C4dbopi2;
  rvec dBOp, dln_BOp_s, dln_BOp_pi, dln_BOp_pi2;
  double *CdboReduction;
} bond_order_data;

typedef struct {
  int nbr;
  int sym_index;
  int dbond_index;
  ivec rel_box;
  //  rvec ext_factor;
  double d;
  rvec dvec;
  bond_order_data bo_data;
} bond_data;


typedef struct {
  int j;
  double val;
} sparse_matrix_entry;

typedef struct {
  int cap, n, m;
  int *start, *end;
  sparse_matrix_entry *entries;
} sparse_matrix;


typedef struct {
  int num_far;
  int H, Htop;
  int hbonds, num_hbonds;
  int bonds, num_bonds;
  int num_3body;
  int gcell_atoms;
} reallocate_data;


typedef struct
{
  int allocated;

  /* communication storage */
  double *tmp_dbl[REAX_MAX_NBRS];
  rvec *tmp_rvec[REAX_MAX_NBRS];
  rvec2 *tmp_rvec2[REAX_MAX_NBRS];
  int  *within_bond_box;

  /* bond order related storage */
  double *total_bond_order;
  double *Deltap, *Deltap_boc;
  double *Delta, *Delta_lp, *Delta_lp_temp, *Delta_e, *Delta_boc, *Delta_val;
  double *dDelta_lp, *dDelta_lp_temp;
  double *nlp, *nlp_temp, *Clp, *vlpex;
  rvec *dDeltap_self;
  int *bond_mark, *done_after;

  /* QEq storage */
  sparse_matrix *H, *L, *U;
  double *Hdia_inv, *b_s, *b_t, *b_prc, *b_prm, *s, *t;
  double *droptol;
  rvec2 *b, *x;

  /* GMRES storage */
  double *y, *z, *g;
  double *hc, *hs;
  double **h, **v;
  /* CG storage */
  double *r, *d, *q, *p;
  rvec2 *r2, *d2, *q2, *p2;
  /* Taper */
  double Tap[8]; //Tap7, Tap6, Tap5, Tap4, Tap3, Tap2, Tap1, Tap0;

  /* storage for analysis */
  int  *mark, *old_mark;
  rvec *x_old;

  /* storage space for bond restrictions */
  int  *restricted;
  int **restricted_list;

  /* integrator */
  rvec *v_const;

  /* force calculations */
  double *CdDelta;  // coefficient of dDelta
  rvec *f;

  /* omp */
  rvec *forceReduction;
  rvec *my_ext_pressReduction;
  double *CdDeltaReduction;
  int *valence_angle_atom_myoffset;

  reallocate_data realloc;
} storage;


typedef union
{
  void *v;
  three_body_interaction_data *three_body_list;
  bond_data          *bond_list;
  dbond_data         *dbo_list;
  dDelta_data        *dDelta_list;
  far_neighbor_data  *far_nbr_list;
  hbond_data         *hbond_list;
} list_type;


struct _reax_list
{
  int allocated;

  int n;
  int num_intrs;

  int *index;
  int *end_index;

  int type;
  list_type select;
  class LAMMPS_NS::Error     *error_ptr;
};
typedef _reax_list  reax_list;


typedef struct
{
  FILE *strj;
  int   trj_offset;
  int   atom_line_len;
  int   bond_line_len;
  int   angle_line_len;
  int   write_atoms;
  int   write_bonds;
  int   write_angles;
  char *line;
  int   buffer_len;
  char *buffer;

  FILE *out;
  FILE *pot;
  FILE *log;
  FILE *mol, *ign;
  FILE *dpl;
  FILE *drft;
  FILE *pdb;
  FILE *prs;

  int   write_steps;
  int   traj_compress;
  int   traj_method;
  char  traj_title[81];
  int   atom_info;
  int   bond_info;
  int   angle_info;

  int   restart_format;
  int   restart_freq;
  int   debug_level;
  int   energy_update_freq;

} output_controls;


typedef struct
{
  int atom_count;
  int atom_list[REAX_MAX_MOLECULE_SIZE];
  int mtypes[REAX_MAX_ATOM_TYPES];
} molecule;


struct LR_data
{
  double H;
  double e_vdW, CEvd;
  double e_ele, CEclmb;

  LAMMPS_INLINE
  LR_data() {}

  LAMMPS_INLINE
  void operator = (const LR_data& rhs) {
    H      = rhs.H;
    e_vdW  = rhs.e_vdW;
    CEvd   = rhs.CEvd;
    e_ele  = rhs.e_ele;
    CEclmb = rhs.CEclmb;
  }
  LAMMPS_INLINE
  void operator = (const LR_data& rhs) volatile {
    H      = rhs.H;
    e_vdW  = rhs.e_vdW;
    CEvd   = rhs.CEvd;
    e_ele  = rhs.e_ele;
    CEclmb = rhs.CEclmb;
  }
};


struct cubic_spline_coef
{
  double a, b, c, d;

  LAMMPS_INLINE
  cubic_spline_coef() {}

  LAMMPS_INLINE
  void operator = (const cubic_spline_coef& rhs) {
    a = rhs.a;
    b = rhs.b;
    c = rhs.c;
    d = rhs.d;
  }
  LAMMPS_INLINE
  void operator = (const cubic_spline_coef& rhs) volatile {
    a = rhs.a;
    b = rhs.b;
    c = rhs.c;
    d = rhs.d;
  }
};



typedef struct
{
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
} LR_lookup_table;
extern LR_lookup_table **LR;

/* function pointer defs */
typedef void (*evolve_function)(reax_system*, control_params*,
                                simulation_data*, storage*, reax_list**,
                                output_controls*, mpi_datatypes* );

typedef void (*interaction_function) (reax_system*, control_params*,
                                      simulation_data*, storage*,
                                      reax_list**, output_controls*);

typedef void (*print_interaction)(reax_system*, control_params*,
                                  simulation_data*, storage*,
                                  reax_list**, output_controls*);

typedef double (*lookup_function)(double);

typedef void (*message_sorter) (reax_system*, int, int, int, mpi_out_data*);
typedef void (*unpacker) ( reax_system*, int, void*, int, neighbor_proc*, int );

typedef void (*dist_packer) (void*, mpi_out_data*);
typedef void (*coll_unpacker) (void*, void*, mpi_out_data*);
#endif

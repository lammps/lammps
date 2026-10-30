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
// clang-format off

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Please cite the related publication:
   H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
   "Parallel Reactive Molecular Dynamics: Numerical Methods and
   Algorithmic Techniques", Parallel Computing, in press.

   QTPIE and QEq with relative electric potentials (fix qtpie/reaxff and
   fix qeq/rel/reaxff): Efstratios M Kritikos (Caltech) and
   Navraj S Lalli (Imperial College London)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(qeq/reaxff,FixQEqReaxFF);
FixStyle(qeq/reax,FixQEqReaxFF);
FixStyle(qtpie/reaxff,FixQEqReaxFF);
FixStyle(qeq/rel/reaxff,FixQEqReaxFF);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_REAXFF_H
#define LMP_FIX_QEQ_REAXFF_H

#include "fix.h"
#include "xl_qeq.h"

namespace LAMMPS_NS {

class FixQEqReaxFF : public Fix {
 public:
  FixQEqReaxFF(class LAMMPS *, int, char **);
  ~FixQEqReaxFF() override;
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  virtual void init_storage();
  void setup_pre_force(int) override;
  void pre_force(int) override;

  void setup_pre_force_respa(int, int) override;
  void pre_force_respa(int, int, int) override;

  void min_setup_pre_force(int);
  void min_pre_force(int) override;

  double compute_scalar() override;

 protected:
  // how the electronegativities in the charge equilibration are determined:
  // per atom type (qeq/reaxff), or as effective electronegativities from
  // overlap-weighted averages (qtpie/reaxff), or with relative potentials of
  // an external electric field (qeq/rel/reaxff)
  enum { QEQ, QTPIE, QEQREL };
  int chi_mode;

  int nevery, reaxflag;
  int matvecs;
  int nn, m_fill;
  int n_cap, nmax, m_cap;
  int pack_flag;
  int nlevels_respa;
  int matrix_free;    // kokkos backend; whether or not to use the matrix-free form
  class NeighList *list;
  class NeighList *list_olap;    // full list for the overlap integrals (QTPIE and QEQREL)
  class PairReaxFF *reaxff;
  class FixEfield *efield;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double swa, swb;     // lower/upper Taper cutoff radius
  double Tap[8];       // Taper function
  double tolerance;    // tolerance for the norm of the rel residual in CG

  double *chi, *eta, *gamma;    // qeq parameters
  double **shld;

  // fictitious charges

  double *s, *t;
  double **s_hist, **t_hist;
  int nprev;

  // extended-Lagrangian (XL) propagation of the fictitious charges

  XLQEq xl;                         // settings and warm-up state
  int xl_active;                    // 1 if the current solve is seeded from the XL history
  double **xls_hist, **xlt_hist;    // auxiliary variable trajectories, newest first
  virtual void xl_update();

  // NOLINTBEGIN
  typedef struct {
    int n, m;
    int *firstnbr;
    int *numnbrs;
    int *jlist;
    double *val;
  } sparse_matrix;
  // NOLINTEND

  sparse_matrix H;
  double *Hdia_inv;
  double *b_s, *b_t;
  double *b_prc, *b_prm;
  double *chi_field;
  double *chi_eff;    // effective electronegativities for the QTPIE and QEQREL modes

  //CG storage
  double *p, *q, *r, *d;
  int imax, maxwarn;

  char *pertype_option;    // argument to determine how per-type info is obtained
  char *gauss_file;        // file with gaussian orbital exponents (QTPIE and QEQREL modes)
  double *gauss_exp;       // gaussian orbital exponents for each atom type
  double **prefactor;      // factor used in computation of overlap integrals
  double **expfactor;      // factor used in exponential term of overlap integrals
  double dist_cutoff_sq;   // squared distance beyond which overlap integrals are neglected
  double scale;            // scaling factor for electric polarization effects

  virtual void pertype_parameters(char *);
  void read_gauss_file();
  void init_olap();
  void init_shielding();
  void init_taper();
  virtual void allocate_storage();
  virtual void deallocate_storage();
  void reallocate_storage();
  virtual void allocate_matrix();
  virtual void deallocate_matrix();
  void reallocate_matrix();

  virtual void init_matvec();
  void init_H();
  virtual void compute_H();
  void compute_H_offsets();
  void compute_H_row(int);
  double calculate_H(double, double);
  virtual void calculate_Q();

  virtual int CG(double *, double *);
  virtual void sparse_matvec(sparse_matrix *, double *, double *);
  virtual void solve_st();

  double rhs_s(int) const;

 public:
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

 protected:
  virtual double parallel_norm(double *, int);
  virtual double parallel_dot(double *, double *, int);
  virtual double parallel_vector_acc(double *, int);

  virtual void vector_sum(double *, double, double *, double, double *, int);
  virtual void vector_add(double *, double, double *, int);

  virtual void get_chi_field();
  virtual void update_chi();
  virtual void calc_chi_eff();
  void prepare_efield_phi();
  void set_efield_phi();
  bool efield_has_variables() const;
  double efield_phi(int) const;
  double chi_eff_one(int) const;

  // settings of the external electric potential for the effective electronegativities
  double phi_ex, phi_ey, phi_ez;
  bool phi_atomstyle;

  // dual CG support
  int dual_enabled;            // 0: Original, separate s & t optimization; 1: dual optimization
  int dual_supported;          // 1 if this variant supports the dual CG optimization
  int matvecs_s, matvecs_t;    // Iteration count for each system
  double cg_resid;             // relative residual at the end of the last CG solve
  double xl_resid;             // larger relative residual of the last s and t solves
};

}    // namespace LAMMPS_NS

#endif
#endif

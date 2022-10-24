/* -*- c++ -*- ---------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_FIX_QEQ_H
#define LMP_FIX_QEQ_H

#include "fix.h"

#define EV_TO_KCAL_PER_MOL 14.4
#define DANGER_ZONE 0.90
#define MIN_CAP 50
#define SAFE_ZONE 1.2
#define MIN_NBRS 100

namespace LAMMPS_NS {

class FixQEq : public Fix {
 public:
  FixQEq(class LAMMPS *, int, char **);
  ~FixQEq() override;
  int setmask() override;
  void init_list(int, class NeighList *) override;
  void setup_pre_force(int) override;
  void setup_pre_force_respa(int, int) override;
  void pre_force_respa(int, int, int) override;
  void min_pre_force(int) override;

  double compute_scalar() override;

  // derived child classes must provide these functions

  void init() override;
  void pre_force(int) override = 0;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  double memory_usage() override;

 protected:
  int nevery;
  int nlocal, nall, m_fill;
  int n_cap, nmax, m_cap;
  int pack_flag;
  int nlevels_respa;
  class NeighList *list;

  int matvecs;
  double qeq_time;

  double swa, swb;             // lower/upper Taper cutoff radius
  double Tap[8];               // Taper function
  double tolerance;            // tolerance for the norm of the rel residual in CG
  int maxiter;                 // maximum number of QEq iterations
  int maxwarn;                 // print warning when max iterations was reached
  double cutoff, cutoff_sq;    // neighbor cutoff

  double *chi, *eta, *gamma, *zeta, *zcore;    // qeq parameters
  double *chizj;
  double **shld;
  int streitz_flag, reax_flag;

  bigint ngroup;

  // fictitious charges

  double *s, *t;
  double **s_hist, **t_hist;
  int nprev;

  typedef struct {
    int n, m;
    int *firstnbr;
    int *numnbrs;
    int *jlist;
    double *val;
  } sparse_matrix;

  sparse_matrix H;
  double *Hdia_inv;
  double *b_s, *b_t;
  double *p, *q, *r, *d;

  // streitz-mintmire

  double alpha;

  // damped dynamics

  double *qf, *q1, *q2, qdamp, qstep;

  // fire
  double *qv;

  void calculate_Q();

  double parallel_norm(double *, int);
  double parallel_dot(double *, double *, int);
  double parallel_vector_acc(double *, int);

  void vector_sum(double *, double, double *, double, double *, int);
  void vector_add(double *, double, double *, int);

  void init_storage();
  void read_file(char *);
  void allocate_storage();
  void deallocate_storage();
  void reallocate_storage();
  void allocate_matrix();
  void deallocate_matrix();
  void reallocate_matrix();

  virtual int CG(double *, double *);
  virtual void sparse_matvec(sparse_matrix *, double *, double *);
};

}    // namespace LAMMPS_NS

#endif

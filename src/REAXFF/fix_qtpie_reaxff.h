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

#ifdef FIX_CLASS
// clang-format off
FixStyle(qtpie/reaxff,FixQtpieReaxFF);
// clang-format on
#else

#ifndef LMP_FIX_QTPIE_REAXFF_H
#define LMP_FIX_QTPIE_REAXFF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixQtpieReaxFF : public Fix {
 public:
  FixQtpieReaxFF(class LAMMPS *, int, char **);
  ~FixQtpieReaxFF() override;
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void init_storage();
  void setup_pre_force(int) override;
  void pre_force(int) override;

  void setup_pre_force_respa(int, int) override;
  void pre_force_respa(int, int, int) override;

  void min_setup_pre_force(int);
  void min_pre_force(int) override;

  double compute_scalar() override;

 protected:
  int nevery, reaxflag;
  int matvecs;
  int nn, nt, m_fill;
  int n_cap, nmax, m_cap;
  int pack_flag;
  int nlevels_respa;
  class NeighList *list;
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
  double *b_prc, *b_prm;
  double *chi_eff;         // array of effective electronegativities

  //CG storage
  double *p, *q, *r, *d;
  int imax, maxwarn;

  char *pertype_option;    // argument to determine how per-type info is obtained
  char *gauss_file;        // input file for gaussian orbital exponents
  double *gauss_exp;       // array of gaussian orbital exponents for each atom type
  double dist_cutoff;      // separation distance beyond which to neglect overlap integrals

  void pertype_parameters(char *);
  void init_shielding();
  void init_taper();
  void allocate_storage();
  void deallocate_storage();
  void reallocate_storage();
  void allocate_matrix();
  void deallocate_matrix();
  void reallocate_matrix();

  void init_matvec();
  void init_H();
  void compute_H();
  double calculate_H(double, double);
  void calculate_Q();

  int CG(double *, double *);
  void sparse_matvec(sparse_matrix *, double *, double *);

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  double parallel_norm(double *, int);
  double parallel_dot(double *, double *, int);
  double parallel_vector_acc(double *, int);

  void vector_sum(double *, double, double *, double, double *, int);
  void vector_add(double *, double, double *, int);

  void calc_chi_eff();
  double find_min_exp(const double*, const int);
  double distance(const double*, const double*);

  int matvecs_s, matvecs_t;    // Iteration count for each system
};

}    // namespace LAMMPS_NS

#endif
#endif

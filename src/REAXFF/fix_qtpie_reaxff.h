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
FixStyle(qtpie/reax,FixQtpieReaxFF);
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
  virtual void init_storage();
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
  int nn, m_fill;
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
  double *chi_field;

  //CG storage
  double *p, *q, *r, *d;
  int imax, maxwarn;

  char *pertype_option;    // argument to determine how per-type info is obtained
			   
  // Params from Kritikos - could rename or move to protected later
  char *gauss_file; // input file for gaussian exponents for each type of REAXFF file
  double cutghost; // ghost atoms cutoff (used for check)
  int nn_prev; // number of local atoms; needed for memory reallocation of chi_eff (when multiprocessing)
  double *gauss_exp; // array of gaussian exponents
  double *chi_eff; // array of effective electronegativities
  double *chi_eff_init; // array of effective electronegativities for FixQEqReax::init_storage()

  // void calculate_chi_eff(LAMMPS_NS::Atom *atom, reax_system *system, double *chi,
  //                        int ni, int nj, double *lchi_eff);
  virtual void pertype_parameters(char *);
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
  double calculate_H(double, double);
  virtual void calculate_Q();

  virtual int CG(double *, double *);
  virtual void sparse_matvec(sparse_matrix *, double *, double *);

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  virtual double parallel_norm(double *, int);
  virtual double parallel_dot(double *, double *, int);
  virtual double parallel_vector_acc(double *, int);

  virtual void vector_sum(double *, double, double *, double, double *, int);
  virtual void vector_add(double *, double, double *, int);

  virtual void get_chi_field();

  // dual CG support
  int dual_enabled;            // 0: Original, separate s & t optimization; 1: dual optimization
  int matvecs_s, matvecs_t;    // Iteration count for each system
};

}    // namespace LAMMPS_NS

#endif
#endif

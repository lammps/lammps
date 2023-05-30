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
FixStyle(acks2/reax,FixACKS2ReaxFF);
FixStyle(acks2/reaxff,FixACKS2ReaxFF);
// clang-format on
#else

#ifndef LMP_FIX_ACKS2_REAXFF_H
#define LMP_FIX_ACKS2_REAXFF_H

#include "fix_qeq_reaxff.h"

namespace LAMMPS_NS {

class FixACKS2ReaxFF : public FixQEqReaxFF {
 public:
  FixACKS2ReaxFF(class LAMMPS *, int, char **);
  ~FixACKS2ReaxFF() override;
  void post_constructor() override;
  void init() override;
  void init_storage() override;
  void pre_force(int) override;

  double *get_s() { return s; }

 protected:
  int NN, last_rows_rank, last_rows_flag;

  double **s_hist_X, **s_hist_last;
  double *bcut_acks2, bond_softness, **bcut;    // acks2 parameters

  sparse_matrix X;
  double *Xdia_inv;
  double *X_diag;

  //BiCGStab storage
  double *g, *q_hat, *r_hat, *y, *z;

  void pertype_parameters(char *) override;
  void init_bondcut();
  void allocate_storage() override;
  void deallocate_storage() override;
  void allocate_matrix() override;
  void deallocate_matrix() override;

  void init_matvec() override;
  void compute_X();
  double calculate_X(double, double);
  void calculate_Q() override;

  int BiCGStab(double *, double *);
  void sparse_matvec_acks2(sparse_matrix *, sparse_matrix *, double *, double *);

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  void more_forward_comm(double *);
  void more_reverse_comm(double *);
  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  double parallel_norm(double *, int) override;
  double parallel_dot(double *, double *, int) override;
  double parallel_vector_acc(double *, int) override;

  void vector_sum(double *, double, double *, double, double *, int) override;
  void vector_add(double *, double, double *, int) override;
  void vector_copy(double *, double *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif

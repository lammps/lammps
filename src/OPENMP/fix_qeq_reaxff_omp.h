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
FixStyle(qeq/reaxff/omp,FixQEqReaxFFOMP);
FixStyle(qeq/reax/omp,FixQEqReaxFFOMP);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_REAXFF_OMP_H
#define LMP_FIX_QEQ_REAXFF_OMP_H

#include "fix_qeq_reaxff.h"

namespace LAMMPS_NS {

class FixQEqReaxFFOMP : public FixQEqReaxFF {

 public:
  FixQEqReaxFFOMP(class LAMMPS *, int, char **);
  ~FixQEqReaxFFOMP() override;
  void init() override;
  void init_storage() override;
  void pre_force(int) override;
  void post_constructor() override;

 protected:
  double **b_temp;

  int do_aspc;
  int aspc_order, aspc_order_max;
  double aspc_omega;
  double *aspc_b;

  void allocate_storage() override;
  void deallocate_storage() override;
  void init_matvec() override;
  void compute_H() override;

  int CG(double *, double *) override;
  void sparse_matvec(sparse_matrix *, double *, double *) override;
  void calculate_Q() override;

  void vector_sum(double *, double, double *, double, double *, int) override;
  void vector_add(double *, double, double *, int) override;

  // dual CG support
  virtual int dual_CG(double *, double *, double *, double *);
  virtual void dual_sparse_matvec(sparse_matrix *, double *, double *, double *);
  virtual void dual_sparse_matvec(sparse_matrix *, double *, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif

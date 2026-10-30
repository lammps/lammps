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
FixStyle(acks2/reaxff/omp,FixACKS2ReaxFFOMP);
FixStyle(acks2/reax/omp,FixACKS2ReaxFFOMP);
// clang-format on
#else

#ifndef LMP_FIX_ACKS2_REAXFF_OMP_H
#define LMP_FIX_ACKS2_REAXFF_OMP_H

#include "fix_acks2_reaxff.h"

namespace LAMMPS_NS {

class FixACKS2ReaxFFOMP : public FixACKS2ReaxFF {
 public:
  FixACKS2ReaxFFOMP(class LAMMPS *, int, char **);
  ~FixACKS2ReaxFFOMP() override;
  double memory_usage() override;

 protected:
  double **b_temp;    // per-thread accumulation buffers
  int nmax_temp;

  void allocate_storage() override;
  void deallocate_storage() override;

  void compute_H() override;
  void compute_X() override;    // NOLINT
  void sparse_matvec_acks2(sparse_matrix *, sparse_matrix *, double *, double *) override;

  double parallel_norm(double *, int) override;
  double parallel_dot(double *, double *, int) override;
  double parallel_vector_acc(double *, int) override;

  void vector_sum(double *, double, double *, double, double *, int) override;
  void vector_add(double *, double, double *, int) override;
  void vector_copy(double *, double *, int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif

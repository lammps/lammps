/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(qeq/reax/omp,FixQEqReaxOMP)

#else

#ifndef LMP_FIX_QEQ_REAX_OMP_H
#define LMP_FIX_QEQ_REAX_OMP_H

#include "fix_qeq_reax.h"

namespace LAMMPS_NS {

class FixQEqReaxOMP : public FixQEqReax {

 public:
  FixQEqReaxOMP(class LAMMPS *, int, char **);
  ~FixQEqReaxOMP();
  virtual void init();
  virtual void init_storage();
  virtual void pre_force(int);
  virtual void post_constructor();

 protected:
  double **b_temp;

  int do_aspc;
  int aspc_order, aspc_order_max;
  double aspc_omega;
  double * aspc_b;

  virtual void allocate_storage();
  virtual void deallocate_storage();
  virtual void init_matvec();
  virtual void compute_H();

  virtual int CG(double*,double*);
  virtual void sparse_matvec(sparse_matrix*,double*,double*);
  virtual void calculate_Q();

  virtual void vector_sum(double*,double,double*,double,double*,int);
  virtual void vector_add(double*, double, double*,int);

  // dual CG support
  virtual int dual_CG(double*,double*,double*,double*);
  virtual void dual_sparse_matvec(sparse_matrix*,double*,double*,double*);
  virtual void dual_sparse_matvec(sparse_matrix*,double*,double*);
};

}

#endif
#endif

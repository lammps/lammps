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

FixStyle(qeq/slater,FixQEqSlater)

#else

#ifndef LMP_FIX_QEQ_SLATER_H
#define LMP_FIX_QEQ_SLATER_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqSlater : public FixQEq {
 public:
  FixQEqSlater(class LAMMPS *, int, char **);
  ~FixQEqSlater() {}
  void init();
  void pre_force(int);

 private:
  void init_matvec();
  int CG(double*,double*);
  void sparse_matvec(sparse_matrix*,double*,double*);
  void compute_H();
  double calculate_H(double, double, double, double, double &);
  double calculate_H_wolf(double, double, double, double, double &);
};
}
#endif
#endif

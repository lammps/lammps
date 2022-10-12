/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(qeq/slater,FixQEqSlater);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_SLATER_H
#define LMP_FIX_QEQ_SLATER_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqSlater : public FixQEq {
 public:
  FixQEqSlater(class LAMMPS *, int, char **);

  void init() override;
  void pre_force(int) override;

 private:
  void init_matvec();
  void sparse_matvec(sparse_matrix *, double *, double *) override;
  void compute_H();
  double calculate_H(double, double, double, double, double &);
  double calculate_H_wolf(double, double, double, double, double &);
  void extract_streitz();

  class PairCoulStreitz *streitz;
};
}    // namespace LAMMPS_NS
#endif
#endif

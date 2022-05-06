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
FixStyle(qeq/shielded,FixQEqShielded);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_SHIELDED_H
#define LMP_FIX_QEQ_SHIELDED_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqShielded : public FixQEq {
 public:
  FixQEqShielded(class LAMMPS *, int, char **);

  void init() override;
  void pre_force(int) override;

 private:
  void extract_reax();
  void init_shielding();
  void init_matvec();
  void compute_H();
  double calculate_H(double, double);
};
}    // namespace LAMMPS_NS
#endif
#endif

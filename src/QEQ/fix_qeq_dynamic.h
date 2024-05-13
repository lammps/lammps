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
FixStyle(qeq/dynamic,FixQEqDynamic);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_DYNAMIC_H
#define LMP_FIX_QEQ_DYNAMIC_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqDynamic : public FixQEq {
 public:
  FixQEqDynamic(class LAMMPS *, int, char **);

  void init() override;
  void pre_force(int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 private:
  double compute_eneg();
};

}    // namespace LAMMPS_NS

#endif
#endif

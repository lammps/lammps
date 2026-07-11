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
FixStyle(qeq/point,FixQEqPoint);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_POINT_H
#define LMP_FIX_QEQ_POINT_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqPoint : public FixQEq {
 public:
  FixQEqPoint(class LAMMPS *, int, char **);

  void init() override;
  void pre_force(int) override;

 private:
  void init_matvec();
  void compute_H();
};
}    // namespace LAMMPS_NS
#endif
#endif

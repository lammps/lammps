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
  ~FixQEqPoint() {}
  void init();
  void pre_force(int);

 private:
  void init_matvec();
  void compute_H();
};
}    // namespace LAMMPS_NS
#endif
#endif

/* ERROR/WARNING messages:

E: Fix qeq/point requires atom attribute q

Self-explanatory.

E: Fix qeq/point group has no atoms

Self-explanatory.

W: H matrix size has been exceeded: m_fill=%d H.m=%d\n

This is the size of the matrix.

E: Fix qeq/point has insufficient QEq matrix size

Occurs when number of neighbor atoms for an atom increased too much
during a run.  Increase SAFE_ZONE and MIN_CAP in fix_qeq.h and
recompile.

*/

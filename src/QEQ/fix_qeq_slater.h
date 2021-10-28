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
  ~FixQEqSlater() {}
  void init();
  void pre_force(int);

 private:
  void init_matvec();
  void sparse_matvec(sparse_matrix *, double *, double *);
  void compute_H();
  double calculate_H(double, double, double, double, double &);
  double calculate_H_wolf(double, double, double, double, double &);
  void extract_streitz();

  class PairCoulStreitz *streitz;
};
}    // namespace LAMMPS_NS
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix qeq/slater requires atom attribute q

Self-explanatory.

E: Fix qeq/slater group has no atoms

Self-explanatory.

E: Invalid param file for fix qeq/slater

Zeta value is 0.0.

E: No pair coul/streitz for fix qeq/slater

These commands must be used together.

E: Fix qeq/slater could not extract params from pair coul/streitz

This should not happen unless pair coul/streitz has been altered.

W: H matrix size has been exceeded: m_fill=%d H.m=%d\n

This is the size of the matrix.

E: Fix qeq/slater has insufficient QEq matrix size

Occurs when number of neighbor atoms for an atom increased too much
during a run.  Increase SAFE_ZONE and MIN_CAP in fix_qeq.h and
recompile.

*/

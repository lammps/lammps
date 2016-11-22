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

FixStyle(qeq/fire,FixQEqFire)

#else

#ifndef LMP_FIX_QEQ_FIRE_H
#define LMP_FIX_QEQ_FIRE_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqFire : public FixQEq {
 public:
  FixQEqFire(class LAMMPS *, int, char **);
  ~FixQEqFire() {}
  void init();
  void pre_force(int);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

 private:
  double compute_eneg();

  class PairComb *comb;
  class PairComb3 *comb3;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix qeq/fire requires atom attribute q

Self-explanatory.

E: Fix qeq/fire group has no atoms

Self-explanatory.

W: Fix qeq/fire tolerance may be too small for damped fires

Self-explanatory.

W: Charges did not converge at step %ld: %lg

Self-explanatory.

*/

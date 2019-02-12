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

FixStyle(qeq/shielded,FixQEqShielded)

#else

#ifndef LMP_FIX_QEQ_SHIELDED_H
#define LMP_FIX_QEQ_SHIELDED_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqShielded : public FixQEq {
 public:
  FixQEqShielded(class LAMMPS *, int, char **);
  ~FixQEqShielded() {}
  void init();
  void pre_force(int);

 private:
  void extract_reax();
  void init_shielding();
  void init_matvec();
  void compute_H();
  double calculate_H(double,double);

};
}
#endif
#endif

/* ERROR/WARNING messages:

E: Fix qeq/shielded requires atom attribute q

Self-explanatory.

E: Fix qeq/shielded group has no atoms

Self-explanatory.

E: Invalid param file for fix qeq/shielded

Invalid value of gamma.

W: Fix qeq has non-zero lower Taper radius cutoff

Absolute value must be <= 0.01.

E: Fix qeq has negative upper Taper radius cutoff

Self-explanatory.

W: Fix qeq has very low Taper radius cutoff

Value should typically be >= 5.0.

W: H matrix size has been exceeded: m_fill=%d H.m=%d\n

This is the size of the matrix.

E: Fix qeq/shielded has insufficient QEq matrix size

Occurs when number of neighbor atoms for an atom increased too much
during a run.  Increase SAFE_ZONE and MIN_CAP in fix_qeq.h and
recompile.

*/

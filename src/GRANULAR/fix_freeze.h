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
FixStyle(freeze,FixFreeze);
// clang-format on
#else

#ifndef LMP_FIX_FREEZE_H
#define LMP_FIX_FREEZE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFreeze : public Fix {
 public:
  FixFreeze(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  virtual void post_force(int);
  void post_force_respa(int, int, int);
  double compute_vector(int);

 protected:
  int force_flag;
  double foriginal[3], foriginal_all[3];
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix freeze requires atom attribute torque

The atom style defined does not have this attribute.

E: More than one fix freeze

Only one of these fixes can be defined, since the granular pair
potentials access it.

*/

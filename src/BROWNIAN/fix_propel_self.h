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
FixStyle(propel/self,FixPropelSelf);
// clang-format on
#else

#ifndef LMP_FIX_PROPEL_SELF_H
#define LMP_FIX_PROPEL_SELF_H

#include "fix.h"
namespace LAMMPS_NS {

class FixPropelSelf : public Fix {
 public:
  FixPropelSelf(class LAMMPS *, int, char **);
  virtual ~FixPropelSelf(){};
  void init();
  void post_force(int);
  void setup(int);
  int setmask();

 private:
  double magnitude;
  double sx, sy, sz;
  int mode;

  void post_force_dipole(int);
  void post_force_velocity(int);
  void post_force_quaternion(int);

  class AtomVecEllipsoid *avec;
};
}    // namespace LAMMPS_NS
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix propel/self command.

Wrong number/type of input arguments.

E: Fix propel/self requires atom attribute mu with option dipole.

Self-explanatory.

E: Fix propel/self requires atom style ellipsoid with option quat.

Self-explanatory.

Fix propel/self requires extended particles with option quat.

Self-explanatory.

*/

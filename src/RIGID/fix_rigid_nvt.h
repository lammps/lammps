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
FixStyle(rigid/nvt,FixRigidNVT);
// clang-format on
#else

#ifndef LMP_FIX_RIGID_NVT_H
#define LMP_FIX_RIGID_NVT_H

#include "fix_rigid_nh.h"

namespace LAMMPS_NS {

class FixRigidNVT : public FixRigidNH {
 public:
  FixRigidNVT(class LAMMPS *, int, char **);
  ~FixRigidNVT() {}
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Did not set temperature for fix rigid/nvt

The temp keyword must be specified.

E: Target temperature for fix rigid/nvt cannot be 0.0

Self-explanatory.

E: Fix rigid/nvt period must be > 0.0

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix rigid/nvt temperature order must be 3 or 5

Self-explanatory.

*/

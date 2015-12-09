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

FixStyle(rigid/npt/small,FixRigidNPTSmall)

#else

#ifndef LMP_FIX_RIGID_NPT_SMALL_H
#define LMP_FIX_RIGID_NPT_SMALL_H

#include "fix_rigid_nh_small.h"

namespace LAMMPS_NS {

class FixRigidNPTSmall : public FixRigidNHSmall {
 public:
  FixRigidNPTSmall(class LAMMPS *, int, char **);
  ~FixRigidNPTSmall() {}
};


}

#endif
#endif

/* ERROR/WARNING messages:

E: Did not set temp or press for fix rigid/npt/small

Self-explanatory.

E: Target temperature for fix rigid/npt/small cannot be 0.0

Self-explanatory.

E: Target pressure for fix rigid/npt/small cannot be < 0.0

Self-explanatory.

E: Fix rigid/npt/small period must be > 0.0

Self-explanatory.

E: Fix rigid npt/small t_chain should not be less than 1

Self-explanatory.

E: Fix rigid npt/small t_order must be 3 or 5

Self-explanatory.

*/

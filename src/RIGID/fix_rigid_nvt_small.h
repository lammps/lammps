/* ----------------------------------------------------------------------
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

FixStyle(rigid/nvt/small,FixRigidNVTSmall)

#else

#ifndef LMP_FIX_RIGID_NVT_SMALL_H
#define LMP_FIX_RIGID_NVT_SMALL_H

#include "fix_rigid_nh_small.h"

namespace LAMMPS_NS {

class FixRigidNVTSmall : public FixRigidNHSmall {
 public:
  FixRigidNVTSmall(class LAMMPS *, int, char **);
  ~FixRigidNVTSmall() {}
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Did not set temp for fix rigid/nvt/small

UNDOCUMENTED

E: Target temperature for fix rigid/nvt/small cannot be 0.0

UNDOCUMENTED

E: Fix rigid/nvt/small period must be > 0.0

UNDOCUMENTED

E: Illegal ... command

UNDOCUMENTED

E: Fix_modify order must be 3 or 5

UNDOCUMENTED

*/

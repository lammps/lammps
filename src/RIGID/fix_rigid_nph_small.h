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

FixStyle(rigid/nph/small,FixRigidNPHSmall)

#else

#ifndef LMP_FIX_RIGID_NPH_SMALL_H
#define LMP_FIX_RIGID_NPH_SMALL_H

#include "fix_rigid_nh_small.h"

namespace LAMMPS_NS {

class FixRigidNPHSmall : public FixRigidNHSmall {
 public:
  FixRigidNPHSmall(class LAMMPS *, int, char **);
  ~FixRigidNPHSmall() {}
};


}

#endif
#endif

/* ERROR/WARNING messages:

E: Pressure control must be used with fix rigid nph/small

UNDOCUMENTED

E: Temperature control must not be used with fix rigid/nph/small

UNDOCUMENTED

E: Target pressure for fix rigid/nph/small cannot be 0.0

UNDOCUMENTED

*/

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

FixStyle(rigid/npt,FixRigidNPT)

#else

#ifndef LMP_FIX_RIGID_NPT_H
#define LMP_FIX_RIGID_NPT_H

#include "fix_rigid_nh.h"

namespace LAMMPS_NS {

class FixRigidNPT : public FixRigidNH {
 public:
  FixRigidNPT(class LAMMPS *, int, char **);
  ~FixRigidNPT() {}
};


}

#endif
#endif

/* ERROR/WARNING messages:

E: Did not set temp or press for fix rigid/npt

UNDOCUMENTED

E: Target temperature for fix rigid/npt cannot be 0.0

UNDOCUMENTED

E: Target pressure for fix rigid/npt cannot be 0.0

UNDOCUMENTED

E: Fix rigid/npt period must be > 0.0

UNDOCUMENTED

E: Illegal ... command

UNDOCUMENTED

E: Fix_modify order must be 3 or 5

UNDOCUMENTED

*/

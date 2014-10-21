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

FixStyle(rigid/nph,FixRigidNPH)

#else

#ifndef LMP_FIX_RIGID_NPH_H
#define LMP_FIX_RIGID_NPH_H

#include "fix_rigid_nh.h"

namespace LAMMPS_NS {

class FixRigidNPH : public FixRigidNH {
 public:
  FixRigidNPH(class LAMMPS *, int, char **);
  ~FixRigidNPH() {}
};


}

#endif
#endif

/* ERROR/WARNING messages:

E: Did not set pressure for fix rigid/nph

The press keyword must be specified.

E: Cannot set temperature for fix rigid/nph

The temp keyword cannot be specified.

*/

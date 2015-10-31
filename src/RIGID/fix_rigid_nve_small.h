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

FixStyle(rigid/nve/small,FixRigidNVESmall)

#else

#ifndef LMP_FIX_RIGID_NVE_SMALL_H
#define LMP_FIX_RIGID_NVE_SMALL_H

#include "fix_rigid_nh_small.h"

namespace LAMMPS_NS {

class FixRigidNVESmall : public FixRigidNHSmall {
 public:
  FixRigidNVESmall(class LAMMPS *, int, char **);
  ~FixRigidNVESmall() {}
};

}

#endif
#endif

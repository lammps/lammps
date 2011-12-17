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

FixStyle(nvt/sphere,FixNVTSphere)

#else

#ifndef LMP_FIX_NVT_SPHERE_H
#define LMP_FIX_NVT_SPHERE_H

#include "fix_nh_sphere.h"

namespace LAMMPS_NS {

class FixNVTSphere : public FixNHSphere {
 public:
  FixNVTSphere(class LAMMPS *, int, char **);
  ~FixNVTSphere() {}
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature control must be used with fix nvt/sphere

Self-explanatory.

E: Pressure control can not be used with fix nvt/sphere

Self-explanatory.

*/

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
FixStyle(nph/sphere/omp,FixNPHSphereOMP);
// clang-format on
#else

#ifndef LMP_FIX_NPH_SPHERE_OMP_H
#define LMP_FIX_NPH_SPHERE_OMP_H

#include "fix_nh_sphere_omp.h"

namespace LAMMPS_NS {

class FixNPHSphereOMP : public FixNHSphereOMP {
 public:
  FixNPHSphereOMP(class LAMMPS *, int, char **);
  virtual ~FixNPHSphereOMP() {}
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature control can not be used with fix nph/sphere

Self-explanatory.

E: Pressure control must be used with fix nph/sphere

Self-explanatory.

*/

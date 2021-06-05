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
FixStyle(nve/dot,FixNVEDot);
// clang-format on
#else

#ifndef LMP_FIX_NVE_DOT_H
#define LMP_FIX_NVE_DOT_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVEDot : public FixNVE {
 public:
  FixNVEDot(class LAMMPS *, int, char **);
  void init();
  void initial_integrate(int);
  void final_integrate();

 private:
  double dt, dthlf, dthlfm;
  class AtomVecEllipsoid *avec;
  // conversion from 3-vector in space frame to 4-vector in body frame
  inline void vec3_to_vec4(const double *q, const double *v3, double *v4)
  {
    v4[0] = -q[1] * v3[0] - q[2] * v3[1] - q[3] * v3[2];
    v4[1] = q[0] * v3[0] + q[3] * v3[1] - q[2] * v3[2];
    v4[2] = -q[3] * v3[0] + q[0] * v3[1] + q[1] * v3[2];
    v4[3] = q[2] * v3[0] - q[1] * v3[1] + q[0] * v3[2];
  }
  // conversion from 4-vector in body frame to 3-vector in space frame
  inline void vec4_to_vec3(const double *q, const double *v4, double *v3)
  {
    v3[0] = -q[1] * v4[0] + q[0] * v4[1] - q[3] * v4[2] + q[2] * v4[3];
    v3[1] = -q[2] * v4[0] + q[3] * v4[1] + q[0] * v4[2] - q[1] * v4[3];
    v3[2] = -q[3] * v4[0] - q[2] * v4[1] + q[1] * v4[2] + q[0] * v4[3];
  }
};

}    // namespace LAMMPS_NS
#endif
#endif

/* ERROR/WARNING messages:

E: Compute nve/dot requires atom style ellipsoid

Self-explanatory.

E: Fix nve/dot requires extended particles

This fix can only be used for particles with a shape setting.

*/

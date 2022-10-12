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
FixStyle(nve/dotc/langevin,FixNVEDotcLangevin);
// clang-format on
#else

#ifndef LMP_FIX_NVE_DOTC_LANGEVIN_H
#define LMP_FIX_NVE_DOTC_LANGEVIN_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVEDotcLangevin : public FixNVE {
 public:
  FixNVEDotcLangevin(class LAMMPS *, int, char **);
  ~FixNVEDotcLangevin() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;

 private:
  double dt, dthlf, dthlfm, dtqrt;
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

 protected:
  int seed;
  class AtomVecEllipsoid *avec;
  double t_start, t_stop, t_period, t_target, tsqrt;
  double gamma, Gamma, ascale;
  double M, gfactor1, gfactor2;
  double gfactor3[3], gfactor4[3], gfactor5[3];
  class RanMars *random;
  void compute_target();
};

}    // namespace LAMMPS_NS
#endif
#endif

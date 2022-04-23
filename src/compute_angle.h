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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(angle,ComputeAngle);
// clang-format on
#else

#ifndef LMP_COMPUTE_ANGLE_H
#define LMP_COMPUTE_ANGLE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeAngle : public Compute {
 public:
  ComputeAngle(class LAMMPS *, int, char **);
  ~ComputeAngle() override;
  void init() override;
  void compute_vector() override;

 private:
  int nsub;
  class AngleHybrid *angle;
  double *emine;
};

}    // namespace LAMMPS_NS

#endif
#endif

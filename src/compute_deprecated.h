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
// list all deprecated and removed compute styles here
ComputeStyle(DEPRECATED,ComputeDeprecated);
// clang-format on
#else

#ifndef LMP_COMPUTE_DEPRECATED_H
#define LMP_COMPUTE_DEPRECATED_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDeprecated : public Compute {
 public:
  ComputeDeprecated(class LAMMPS *, int, char **);
  void init() override {}
};

}    // namespace LAMMPS_NS

#endif
#endif

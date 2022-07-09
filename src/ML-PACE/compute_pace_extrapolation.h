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

//
// Created by Yury Lysogorskiy on 23.06.22.
//

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(pace/extrapolation,ComputePACEExtrapolation);
// clang-format on
#else

#ifndef COMPUTE_PACE_H
#define COMPUTE_PACE_H

#include "compute.h"
#include "pair_pace_extrapolation.h"

namespace LAMMPS_NS {
class PairPACEExtrapolation;

class ComputePACEExtrapolation : public Compute {
 public:
  ComputePACEExtrapolation(class LAMMPS *, int, char **);
  void init() override;
  double compute_scalar() override;
  void compute_peratom() override;

 private:
  PairPACEExtrapolation *pair_pace_extrapolation;
  void invoke_compute_extrapolation_grades();
};

}    // namespace LAMMPS_NS
#endif    //COMPUTE_PACE_H
#endif

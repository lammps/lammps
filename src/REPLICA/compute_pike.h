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
ComputeStyle(pike,ComputePIKE);
// clang-format on
#else

#ifndef LMP_COMPUTE_PIKE_H
#define LMP_COMPUTE_PIKE_H

#include "compute.h"
#include "fix_dp_pimd.h"

namespace LAMMPS_NS {

class ComputePIKE : public Compute {
 public:
  ComputePIKE(class LAMMPS *, int, char **);
  void init();
  double compute_scalar();

 private:
  double pfactor;
  int bead_flag;
  char *id_fix;
  int ifix;
  class FixDPPimd *fix_dppimd;
};

}    // namespace LAMMPS_NS

#endif
#endif
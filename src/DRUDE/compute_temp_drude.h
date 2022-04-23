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
ComputeStyle(temp/drude,ComputeTempDrude);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_DRUDE_H
#define LMP_COMPUTE_TEMP_DRUDE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempDrude : public Compute {
 public:
  ComputeTempDrude(class LAMMPS *, int, char **);
  ~ComputeTempDrude() override;
  void init() override;
  void setup() override;
  void compute_vector() override;
  double compute_scalar() override;
  int modify_param(int, char **);

 private:
  class FixDrude *fix_drude;
  char *id_temp;
  class Compute *temperature;
  bigint dof_core, dof_drude;
  double kineng_core, kineng_drude;
  double temp_core, temp_drude;

  void dof_compute();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(temp/ramp,ComputeTempRamp);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_RAMP_H
#define LMP_COMPUTE_TEMP_RAMP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempRamp : public Compute {
 public:
  ComputeTempRamp(class LAMMPS *, int, char **);
  ~ComputeTempRamp() override;
  void init() override {}
  void setup() override;
  double compute_scalar() override;
  void compute_vector() override;

  void remove_bias(int, double *) override;
  void remove_bias_all() override;
  void remove_bias_thr(int, double *, double *) override;
  void restore_bias(int, double *) override;
  void restore_bias_thr(int, double *, double *) override;
  void restore_bias_all() override;
  double memory_usage() override;

 private:
  int coord_dim;
  double coord_lo, coord_hi;
  int v_dim;
  double v_lo, v_hi;
  int scaleflag;
  double tfactor, xscale, yscale, zscale;

  void dof_compute();
};

}    // namespace LAMMPS_NS

#endif
#endif

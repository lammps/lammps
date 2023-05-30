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
ComputeStyle(temp/partial,ComputeTempPartial);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_PARTIAL_H
#define LMP_COMPUTE_TEMP_PARTIAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempPartial : public Compute {
 public:
  ComputeTempPartial(class LAMMPS *, int, char **);
  ~ComputeTempPartial() override;
  void init() override {}
  void setup() override;
  double compute_scalar() override;
  void compute_vector() override;

  int dof_remove(int) override;
  void remove_bias(int, double *) override;
  void remove_bias_thr(int, double *, double *) override;
  void remove_bias_all() override;
  void reapply_bias_all() override;
  void restore_bias(int, double *) override;
  void restore_bias_thr(int, double *, double *) override;
  void restore_bias_all() override;
  double memory_usage() override;

 protected:
  int xflag, yflag, zflag;
  double tfactor;

  void dof_compute();
};

}    // namespace LAMMPS_NS

#endif
#endif

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
ComputeStyle(temp/sphere,ComputeTempSphere);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_SPHERE_H
#define LMP_COMPUTE_TEMP_SPHERE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempSphere : public Compute {
 public:
  ComputeTempSphere(class LAMMPS *, int, char **);
  ~ComputeTempSphere() override;
  void init() override;
  void setup() override;
  double compute_scalar() override;
  void compute_vector() override;

  void remove_bias(int, double *) override;
  void remove_bias_thr(int, double *, double *) override;
  void restore_bias(int, double *) override;
  void restore_bias_thr(int, double *, double *) override;

 private:
  int mode;
  double tfactor;
  char *id_bias;
  Compute *tbias;    // ptr to additional bias compute

  void dof_compute();
};

}    // namespace LAMMPS_NS

#endif
#endif

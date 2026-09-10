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
ComputeStyle(temp/deform,ComputeTempDeform);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_DEFORM_H
#define LMP_COMPUTE_TEMP_DEFORM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempDeform : public Compute {
 public:
  ComputeTempDeform(class LAMMPS *, int, char **);
  ~ComputeTempDeform() override;
  void post_constructor() override;
  void init() override;
  void setup() override;
  double compute_scalar() override;
  void compute_vector() override;

  void remove_bias(int, double *) override;
  void remove_bias_thr(int, double *, double *) override;
  void remove_bias_all() override;
  void restore_bias(int, double *) override;
  void restore_bias_thr(int, double *, double *) override;
  void restore_bias_all() override;
  double memory_usage() override;
  int modify_param(int, char **) override;

  void remove_deform_bias(int, double *);
  void remove_deform_bias_thr(int, double *, double *);
  virtual void remove_deform_bias_all();
  virtual void remove_deform_bias_all_kk() {}
  void restore_deform_bias(int, double *);
  void restore_deform_bias_thr(int, double *, double *);
  virtual void restore_deform_bias_all();
  virtual void restore_deform_bias_all_kk() {}
  void apply_deform_bias(double *, double *, double *, double *, double, double);
  virtual void apply_deform_bias_all(double dtv = 0.0);
  virtual void apply_deform_bias_all_kk(double) {}

  class Compute* temperature; // internal temperature compute
  int which;                  // whether internal temp compute has a bias

 protected:
  char *id_temp;
  int tcomputeflag; // 1 = internal temp compute created by this compute

  virtual void dof_compute();
};

}    // namespace LAMMPS_NS

#endif
#endif

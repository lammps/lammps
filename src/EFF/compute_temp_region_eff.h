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
ComputeStyle(temp/region/eff,ComputeTempRegionEff);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_REGION_EFF_H
#define LMP_COMPUTE_TEMP_REGION_EFF_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempRegionEff : public Compute {
 public:
  ComputeTempRegionEff(class LAMMPS *, int, char **);
  virtual ~ComputeTempRegionEff();
  void init();
  void setup();
  virtual double compute_scalar();
  virtual void compute_vector();

  void dof_remove_pre(void);
  int dof_remove(int);
  void remove_bias(int, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_all();
  double memory_usage();

 protected:
  int iregion;
  char *idregion;
};

}    // namespace LAMMPS_NS

#endif
#endif

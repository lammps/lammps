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
ComputeStyle(temp/deform/eff,ComputeTempDeformEff);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_DEFORM_EFF_H
#define LMP_COMPUTE_TEMP_DEFORM_EFF_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempDeformEff : public Compute {
 public:
  ComputeTempDeformEff(class LAMMPS *, int, char **);
  ~ComputeTempDeformEff() override;
  void init() override;
  void setup() override;
  double compute_scalar() override;
  void compute_vector() override;

  void remove_bias(int, double *) override;
  void remove_bias_all() override;
  void restore_bias(int, double *) override;
  void restore_bias_all() override;
  double memory_usage() override;

 protected:
  double tfactor;
  double vbias[3];      // stored velocity bias for one atom
  double **vbiasall;    // stored velocity bias for all atoms
  int maxbias;          // size of vbiasall array

  virtual void dof_compute();
};

}    // namespace LAMMPS_NS

#endif
#endif

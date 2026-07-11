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

#ifdef FIX_CLASS
// clang-format off
FixStyle(rheo/viscosity,FixRHEOViscosity)
// clang-format on
#else

#ifndef LMP_FIX_RHEO_VISCOSITY_H
#define LMP_FIX_RHEO_VISCOSITY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRHEOViscosity : public Fix {
 public:
  FixRHEOViscosity(class LAMMPS *, int, char **);
  ~FixRHEOViscosity() override;
  int setmask() override;
  void init() override;
  void setup_pre_force(int) override;
  void post_neighbor() override;
  void pre_force(int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  double *eta, *npow, *K, *gd0, *tau0;
  int *viscosity_style, constant_flag, evolve_flag;

  class FixRHEO *fix_rheo;
  class ComputeRHEOGrad *compute_grad;
};

}    // namespace LAMMPS_NS

#endif
#endif

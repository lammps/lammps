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

#ifdef FIX_CLASS
// clang-format off
FixStyle(viscous,FixViscous);
// clang-format on
#else

#ifndef LMP_FIX_VISCOUS_H
#define LMP_FIX_VISCOUS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixViscous : public Fix {
 public:
  FixViscous(class LAMMPS *, int, char **);
  ~FixViscous() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;

 protected:
  double *gamma;
  int ilevel_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif

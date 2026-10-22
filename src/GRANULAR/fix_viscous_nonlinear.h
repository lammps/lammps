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
FixStyle(viscous/nonlinear,FixViscousNonlinear);
// clang-format on
#else

#ifndef LMP_FIX_VISCOUS_NONLINEAR_H
#define LMP_FIX_VISCOUS_NONLINEAR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixViscousNonlinear : public Fix {
 public:
  FixViscousNonlinear(class LAMMPS *, int, char **);
  ~FixViscousNonlinear() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;

 protected:
  double rho_fluid, mu_fluid, v_fluid[3];
  int ilevel_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif

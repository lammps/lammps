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
FixStyle(langevin/eff,FixLangevinEff);
// clang-format on
#else

#ifndef LMP_FIX_LANGEVIN_EFF_H
#define LMP_FIX_LANGEVIN_EFF_H

#include "fix_langevin.h"

namespace LAMMPS_NS {

class FixLangevinEff : public FixLangevin {
 public:
  FixLangevinEff(class LAMMPS *, int, char **);
  ~FixLangevinEff() override;
  void end_of_step() override;
  double compute_scalar() override;
  void post_force(int) override;

 private:
  double *erforcelangevin;

  void post_force_no_tally();
  void post_force_tally();
};

}    // namespace LAMMPS_NS

#endif
#endif

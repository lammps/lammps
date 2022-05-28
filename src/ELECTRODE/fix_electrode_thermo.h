/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

// clang-format off
FixStyle(electrode/thermo, FixElectrodeThermo);
// clang-format on

#else

#ifndef LMP_FIX_ELECTRODE_THERMO_H
#define LMP_FIX_ELECTRODE_THERMO_H

#include "fix_electrode_conp.h"

namespace LAMMPS_NS {

class FixElectrodeThermo : public FixElectrodeConp {
 public:
  FixElectrodeThermo(class LAMMPS *, int, char **);
  ~FixElectrodeThermo() override;
  void update_psi() override;
  void pre_update() override;

 protected:
 private:
  void compute_macro_matrices() override;
  class RanMars *thermo_random;
  double delta_psi_0;
  double group_q_old[2];
  double vac_cap;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (GU), Robert Meissner (Hereon, TUHH)
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
  virtual void update_psi_set_constraint() override;

 private:
  class RanMars *thermo_random;
  double delta_v_0; // target voltage
};

}    // namespace LAMMPS_NS

#endif
#endif

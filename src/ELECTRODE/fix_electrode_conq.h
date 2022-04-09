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
FixStyle(electrode/conq, FixElectrodeConq);
// clang-format on

#else

#ifndef LMP_FIX_ELECTRODE_CONQ_H
#define LMP_FIX_ELECTRODE_CONQ_H

#include "fix_electrode_conp.h"

namespace LAMMPS_NS {

class FixElectrodeConq : public FixElectrodeConp {
 public:
  FixElectrodeConq(class LAMMPS *, int, char **);
  void update_psi() override;

 private:
  std::vector<double> group_q;
};

}    // namespace LAMMPS_NS

#endif
#endif

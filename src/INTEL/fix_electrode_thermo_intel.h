/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Shern Tee (UQ)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

// clang-format off
FixStyle(electrode/thermo/intel, FixElectrodeThermoIntel)
// clang-format on

#else

#ifndef LMP_FIX_ELECTRODE_THERMO_INTEL_H
#define LMP_FIX_ELECTRODE_THERMO_INTEL_H

#include "electrode_accel_intel.h"
#include "fix_electrode_thermo.h"

namespace LAMMPS_NS {

class FixElectrodeThermoIntel : public FixElectrodeThermo {
 public:
  FixElectrodeThermoIntel(class LAMMPS *lmp, int narg, char **arg) :
      FixElectrodeThermo(lmp, narg, arg)
  {
    intelflag = true;
    delete accel_interface;
    accel_interface = new ElectrodeAccelIntel(lmp);
  }
};

}    // namespace LAMMPS_NS

#endif
#endif

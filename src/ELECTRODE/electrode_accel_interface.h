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

#ifndef LMP_ELECTRODE_ACCEL_INTERFACE_H
#define LMP_ELECTRODE_ACCEL_INTERFACE_H

#include "pointers.h"

namespace LAMMPS_NS {

class ElectrodeAccelInterface : protected Pointers {
 public:
  ElectrodeAccelInterface(class LAMMPS *lmp) : Pointers(lmp) {}
  virtual void intel_find_fix() {}
  virtual void intel_pack_buffers() {}
};

}    // namespace LAMMPS_NS
#endif

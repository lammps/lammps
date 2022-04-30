/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_ELECTRODE_ACCEL_INTEL_H
#define LMP_ELECTRODE_ACCEL_INTEL_H

#include "electrode_accel_interface.h"
#include "fix_intel.h"
#include "intel_buffers.h"

namespace LAMMPS_NS {

class ElectrodeAccelIntel : public ElectrodeAccelInterface {
 public:
  ElectrodeAccelIntel(class LAMMPS *lmp);
  void intel_find_fix();
  void intel_pack_buffers();

 private:
  class FixIntel *fix;
  template <class flt_t, class acc_t>
  void intel_pack_buffers_prec(IntelBuffers<flt_t, acc_t> *buffers);
};

}    // namespace LAMMPS_NS

#endif

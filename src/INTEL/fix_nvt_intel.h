// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nvt/intel,FixNVTIntel);
// clang-format on
#else

#ifndef LMP_FIX_NVT_INTEL_H
#define LMP_FIX_NVT_INTEL_H

#include "fix_nh_intel.h"

namespace LAMMPS_NS {

class FixNVTIntel : public FixNHIntel {
 public:
  FixNVTIntel(class LAMMPS *, int, char **);
  ~FixNVTIntel() {}
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature control must be used with fix nvt

Self-explanatory.

E: Pressure control can not be used with fix nvt

Self-explanatory.

*/

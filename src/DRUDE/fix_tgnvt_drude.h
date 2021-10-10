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
FixStyle(tgnvt/drude,FixTGNVTDrude);
// clang-format on
#else

#ifndef LMP_FIX_TGNVT_DRUDE_H
#define LMP_FIX_TGNVT_DRUDE_H

#include "fix_tgnh_drude.h"

namespace LAMMPS_NS {

class FixTGNVTDrude : public FixTGNHDrude {
 public:
  FixTGNVTDrude(class LAMMPS *, int, char **);
  ~FixTGNVTDrude() {}
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

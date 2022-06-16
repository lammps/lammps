/* *- c++ -*- -----------------------------------------------------------
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
   Contributing authors: William McDoniel (RWTH Aachen University)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/long/coul/long/intel,PairLJLongCoulLongIntel);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_LONG_COUL_LONG_INTEL_H
#define LMP_PAIR_LJ_LONG_COUL_LONG_INTEL_H

#include "fix_intel.h"
#include "pair_lj_long_coul_long.h"

namespace LAMMPS_NS {
class PairLJLongCoulLongIntel : public PairLJLongCoulLong {
 public:
  PairLJLongCoulLongIntel(class LAMMPS *);

  void init_style() override;
};
}    // namespace LAMMPS_NS
#endif
#endif

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
   Contributing author: Markus Hohnerbach (RWTH)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(airebo/morse/intel,PairAIREBOMorseIntel);
// clang-format on
#else

#ifndef LMP_PAIR_AIREBO_MORSE_INTEL_H
#define LMP_PAIR_AIREBO_MORSE_INTEL_H

#include "pair_airebo_intel.h"

namespace LAMMPS_NS {

class PairAIREBOMorseIntel : public PairAIREBOIntel {
 public:
  PairAIREBOMorseIntel(class LAMMPS *);
  virtual void settings(int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif

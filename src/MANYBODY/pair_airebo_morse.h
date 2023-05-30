/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(airebo/morse,PairAIREBOMorse);
// clang-format on
#else

#ifndef LMP_PAIR_AIREBO_MORSE_H
#define LMP_PAIR_AIREBO_MORSE_H

#include "pair_airebo.h"

namespace LAMMPS_NS {

class PairAIREBOMorse : public PairAIREBO {
 public:
  PairAIREBOMorse(class LAMMPS *);
  void settings(int, char **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif

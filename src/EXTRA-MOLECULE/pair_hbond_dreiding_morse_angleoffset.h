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
PairStyle(hbond/dreiding/morse/angleoffset,PairHbondDreidingMorseAngleoffset);
// clang-format on
#else

#ifndef LMP_PAIR_HBOND_DREIDING_MORSE_ANGLEOFFSET_H
#define LMP_PAIR_HBOND_DREIDING_MORSE_ANGLEOFFSET_H

#include "pair_hbond_dreiding_lj_angleoffset.h"

namespace LAMMPS_NS {

class PairHbondDreidingMorseAngleoffset : public PairHbondDreidingLJangleoffset {
 public:
  PairHbondDreidingMorseAngleoffset(class LAMMPS *);

  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double single(int, int, int, int, double, double, double, double &) override;
};

}    // namespace LAMMPS_NS

#endif
#endif

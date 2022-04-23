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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lubricateU/poly,PairLubricateUPoly);
// clang-format on
#else

#ifndef LMP_PAIR_LUBRICATEU_POLY_H
#define LMP_PAIR_LUBRICATEU_POLY_H

#include "pair_lubricateU.h"

namespace LAMMPS_NS {

class PairLubricateUPoly : public PairLubricateU {
 public:
  PairLubricateUPoly(class LAMMPS *);

  void compute(int, int) override;
  void settings(int, char **) override;
  void init_style() override;

 private:
  double vol_P;
  int flagdeform, flagwall, flagVF, flagHI;
  class FixWall *wallfix;

  void iterate(double **, int);
  void compute_RE(double **) override;
  void compute_RU(double **) override;
  void compute_Fh(double **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif

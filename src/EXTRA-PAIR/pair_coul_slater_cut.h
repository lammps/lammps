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
PairStyle(coul/slater/cut,PairCoulSlaterCut);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_SLATER_CUT_H
#define LMP_PAIR_COUL_SLATER_CUT_H

#include "pair_coul_cut.h"

namespace LAMMPS_NS {

class PairCoulSlaterCut : public PairCoulCut {
 public:
  PairCoulSlaterCut(class LAMMPS *);
  void compute(int, int) override;
  void settings(int, char **) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  double lamda;
};

}    // namespace LAMMPS_NS

#endif
#endif

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
PairStyle(tersoff/mod/c,PairTersoffMODC);
// clang-format on
#else

#ifndef LMP_PAIR_TERSOFF_MOD_C_H
#define LMP_PAIR_TERSOFF_MOD_C_H

#include "pair_tersoff_mod.h"

namespace LAMMPS_NS {

class PairTersoffMODC : public PairTersoffMOD {
 public:
  PairTersoffMODC(class LAMMPS *lmp) : PairTersoffMOD(lmp){};

  static constexpr int NPARAMS_PER_LINE = 21;

 protected:
  void read_file(char *) override;
  void repulsive(Param *, double, double &, int, double &) override;
};

}    // namespace LAMMPS_NS

#endif
#endif

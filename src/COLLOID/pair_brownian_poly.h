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
PairStyle(brownian/poly,PairBrownianPoly);
// clang-format on
#else

#ifndef LMP_PAIR_BROWNIAN_POLY_H
#define LMP_PAIR_BROWNIAN_POLY_H

#include "pair_brownian.h"

namespace LAMMPS_NS {

class PairBrownianPoly : public PairBrownian {
 public:
  PairBrownianPoly(class LAMMPS *);
  ~PairBrownianPoly() {}
  void compute(int, int);
  double init_one(int, int);
  void init_style();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Pair brownian/poly requires newton pair off

Self-explanatory.

E: Pair brownian/poly requires atom style sphere

Self-explanatory.

E: Pair brownian/poly requires extended particles

One of the particles has radius 0.0.

E: Cannot use multiple fix wall commands with pair brownian

Self-explanatory.

*/

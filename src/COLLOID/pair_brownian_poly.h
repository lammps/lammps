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

  void compute(int, int) override;
  double init_one(int, int) override;
  void init_style() override;

 protected:
  // deterministic, order- and MPI-rank-independent uniform random number for
  // a pair, keyed on the unordered atom-tag pair, timestep, seed, and stream
  // index.  Both halves of a pair draw the identical value so the pairwise
  // Brownian force obeys Newton's 3rd law under "newton off".  See issue #2933.
  static double pair_uniform(tagint, tagint, bigint, int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif

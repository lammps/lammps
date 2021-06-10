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
PairStyle(yukawa/colloid,PairYukawaColloid);
// clang-format on
#else

#ifndef LMP_PAIR_YUKAWA_COLLOID_H
#define LMP_PAIR_YUKAWA_COLLOID_H

#include "pair_yukawa.h"

namespace LAMMPS_NS {

class PairYukawaColloid : public PairYukawa {
 public:
  PairYukawaColloid(class LAMMPS *);
  virtual ~PairYukawaColloid() {}
  virtual void compute(int, int);
  void init_style();
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Pair yukawa/colloid requires atom style sphere

Self-explanatory.

E: Pair yukawa/colloid requires atoms with same type have same radius

Self-explanatory.

*/

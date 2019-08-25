/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(DEPRECATED,PairDeprecated)
PairStyle(reax,PairDeprecated)

#else

#ifndef LMP_PAIR_DEPRECATED_H
#define LMP_PAIR_DEPRECATED_H

#include "pair.h"

namespace LAMMPS_NS {

class PairDeprecated : public Pair {
 public:
  PairDeprecated(class LAMMPS *lmp) : Pair(lmp) {}
  virtual ~PairDeprecated() {}

  virtual void compute(int, int) {}
  virtual void settings(int, char **);
  virtual void coeff(int, char **) {}

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

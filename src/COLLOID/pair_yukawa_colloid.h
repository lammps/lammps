/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PAIR_YUKAWA_COLLOID_H
#define PAIR_YUKAWA_COLLOID_H

#include "pair_yukawa.h"

namespace LAMMPS_NS {

class PairYukawaColloid : public PairYukawa {
 public:
  PairYukawaColloid(class LAMMPS *);
  ~PairYukawaColloid() {}
  void compute(int, int);
  void init_style();
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double &);
};

}

#endif

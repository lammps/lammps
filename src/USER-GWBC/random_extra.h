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

#ifndef LMP_RANEXTRA_H
#define LMP_RANEXTRA_H

#include "pointers.h"

namespace LAMMPS_NS {

class RanExtra : protected Pointers {
 public:
  RanExtra(class LAMMPS *, int);
  double uniform();
  double gaussian(double mu, double sigma);
  double rayleigh(double sigma);
  double besselexp(double theta, double alpha, double cp);

  void reset(int);
  int state();

 private:
  int seed,save;
  double second;
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid seed for Park random # generator

The initial seed for this random number generator must be a positive
integer.

*/

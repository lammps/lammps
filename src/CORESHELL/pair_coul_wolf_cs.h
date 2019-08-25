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

PairStyle(coul/wolf/cs,PairCoulWolfCS)

#else

#ifndef LMP_PAIR_COUL_WOLF_CS_H_
#define LMP_PAIR_COUL_WOLF_CS_H_

#include "pair_coul_wolf.h"

namespace LAMMPS_NS {

class PairCoulWolfCS : public PairCoulWolf {
 public:
  PairCoulWolfCS( class LAMMPS *);
  virtual void compute( int, int);
};

}

#endif
#endif /* LMP_PAIR_COUL_WOLF_CS_H_ */

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair coul/wolf/cs requires atom attribute q

The atom style defined does not have this attribute.

*/

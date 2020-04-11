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

PairStyle(lj/class2/coul/long/cs,PairLJClass2CoulLongCS)

#else

#ifndef LMP_PAIR_LJ_CLASS2_COUL_LONG_CS_H
#define LMP_PAIR_LJ_CLASS2_COUL_LONG_CS_H

#include "pair_lj_class2_coul_long.h"

namespace LAMMPS_NS {

class PairLJClass2CoulLongCS : public PairLJClass2CoulLong {

 public:
  PairLJClass2CoulLongCS(class LAMMPS *);
  virtual void compute(int, int);
  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style lj/class2/coul/long requires atom attribute q

The atom style defined does not have this attribute.

E: Pair style requires a KSpace style

No kspace style is defined.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/

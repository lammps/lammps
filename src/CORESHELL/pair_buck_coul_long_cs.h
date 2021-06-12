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
PairStyle(buck/coul/long/cs,PairBuckCoulLongCS);
// clang-format on
#else

#ifndef LMP_PAIR_BUCK_COUL_LONG_CS_H
#define LMP_PAIR_BUCK_COUL_LONG_CS_H

#include "pair_buck_coul_long.h"

namespace LAMMPS_NS {

class PairBuckCoulLongCS : public PairBuckCoulLong {
 public:
  PairBuckCoulLongCS(class LAMMPS *);
  virtual void compute(int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Pair style buck/coul/long requires atom attribute q

The atom style defined does not have these attributes.

E: Pair style requires a KSpace style

No kspace style is defined.

*/

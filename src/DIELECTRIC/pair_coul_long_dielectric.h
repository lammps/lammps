/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(coul/long/dielectric,PairCoulLongDielectric);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_LONG_DIELECTRIC_H
#define LMP_PAIR_COUL_LONG_DIELECTRIC_H

#include "pair_coul_long.h"

namespace LAMMPS_NS {

class PairCoulLongDielectric : public PairCoulLong {
 public:
  PairCoulLongDielectric(class LAMMPS *);
  ~PairCoulLongDielectric();
  virtual void compute(int, int);
  virtual void init_style();

  double **efield;

 protected:
  class AtomVecDielectric *avec;
  int nmax;
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

E: Pair style coul/long/dielectric requires atom attribute q

The atom style defined does not have this attribute.

E: Pair style requires a KSpace style

No kspace style is defined.

*/

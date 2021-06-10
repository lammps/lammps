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
PairStyle(coul/cut/global,PairCoulCutGlobal);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_CUT_GLOBAL_H
#define LMP_PAIR_COUL_CUT_GLOBAL_H

#include "pair_coul_cut.h"

namespace LAMMPS_NS {

class PairCoulCutGlobal : public PairCoulCut {
 public:
  PairCoulCutGlobal(class LAMMPS *lmp) : PairCoulCut(lmp) {}
  void coeff(int, char **);
  void *extract(const char *, int &);
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

E: Pair style coul/cut requires atom attribute q

The atom style defined does not have these attributes.

*/

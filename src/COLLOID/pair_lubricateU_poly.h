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
PairStyle(lubricateU/poly,PairLubricateUPoly);
// clang-format on
#else

#ifndef LMP_PAIR_LUBRICATEU_POLY_H
#define LMP_PAIR_LUBRICATEU_POLY_H

#include "pair_lubricateU.h"

namespace LAMMPS_NS {

class PairLubricateUPoly : public PairLubricateU {
 public:
  PairLubricateUPoly(class LAMMPS *);
  ~PairLubricateUPoly() {}
  void compute(int, int);
  void settings(int, char **);
  void init_style();

 private:
  double vol_P;
  int flagdeform, flagwall, flagVF, flagHI;
  class FixWall *wallfix;

  void iterate(double **, int);
  void compute_RE(double **);
  void compute_RU(double **);
  void compute_Fh(double **);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: Cannot include log terms without 1/r terms; setting flagHI to 1

Self-explanatory.

E: Pair lubricateU/poly requires newton pair off

Self-explanatory.

E: Pair lubricateU/poly requires ghost atoms store velocity

Use the comm_modify vel yes command to enable this.

E: Pair lubricate/poly requires atom style sphere

Self-explanatory.

E: Pair lubricate/poly requires extended particles

One of the particles has radius 0.0.

E: Cannot use multiple fix wall commands with pair lubricateU

Self-explanatory.

*/

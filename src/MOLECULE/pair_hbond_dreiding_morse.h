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

#ifdef PAIR_CLASS

PairStyle(hbond/dreiding/morse,PairHbondDreidingMorse)

#else

#ifndef LMP_PAIR_HBOND_DREIDING_MORSE_H
#define LMP_PAIR_HBOND_DREIDING_MORSE_H

#include "pair_hbond_dreiding_lj.h"

namespace LAMMPS_NS {

class PairHbondDreidingMorse : public PairHbondDreidingLJ {
 public:
  PairHbondDreidingMorse(class LAMMPS *);
  virtual ~PairHbondDreidingMorse() {};
  virtual void compute(int, int);
  void coeff(int, char **);
  void init_style();
  double single(int, int, int, int, double, double, double, double &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair inner cutoff >= Pair outer cutoff

The specified cutoffs for the pair style are inconsistent.

E: Pair style hbond/dreiding requires molecular system

Self-explanatory.

E: Pair style hbond/dreiding requires atom IDs

Self-explanatory.

E: Pair style hbond/dreiding requires an atom map, see atom_modify

Self-explanatory.

E: Pair style hbond/dreiding requires newton pair on

See the newton command for details.

E: No pair hbond/dreiding coefficients set

Self-explanatory.

*/

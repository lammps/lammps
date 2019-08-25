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

PairStyle(born/coul/dsf/cs,PairBornCoulDSFCS)

#else

#ifndef LMP_PAIR_BORN_COUL_DSF_CS_H
#define LMP_PAIR_BORN_COUL_DSF_CS_H

#include "pair_born_coul_dsf.h"

namespace LAMMPS_NS {

class PairBornCoulDSFCS : public PairBornCoulDSF {
 public:
  PairBornCoulDSFCS(class LAMMPS *);
  virtual void compute(int, int);
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

E: Pair style born/coul/dsf requires atom attribute q

The atom style defined does not have this attribute.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/

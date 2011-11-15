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

PairStyle(lubricate/poly,PairLubricatePoly)

#else

#ifndef LMP_PAIR_LUBRICATE_POLY_H
#define LMP_PAIR_LUBRICATE_POLY_H

#include "pair_lubricate.h"

namespace LAMMPS_NS {

class PairLubricatePoly : public PairLubricate {
 public:
  PairLubricatePoly(class LAMMPS *);
  ~PairLubricatePoly() {}
  void compute(int, int);
  void init_style();
};

}

#endif
#endif

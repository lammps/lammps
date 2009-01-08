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

#ifndef PAIR_GRAN_HOOKE_H
#define PAIR_GRAN_HOOKE_H

#include "pair_gran_hooke_history.h"

namespace LAMMPS_NS {

class PairGranHooke : public PairGranHookeHistory {
 public:
  PairGranHooke(class LAMMPS *);
  void compute(int, int);
};

}

#endif

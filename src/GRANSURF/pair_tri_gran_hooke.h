/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(tri/gran/hooke,PairTriGranHooke)

#else

#ifndef LMP_PAIR_TRI_GRAN_HOOKE_H
#define LMP_PAIR_TRI_GRAN_HOOKE_H

#include "pair_tri_gran_hooke_history.h"

namespace LAMMPS_NS {

class PairTriGranHooke : public PairTriGranHookeHistory {
 public:
  PairTriGranHooke(class LAMMPS *);
  ~PairTriGranHooke() {}
  void compute(int, int);
};

}

#endif
#endif

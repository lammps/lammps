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
// clang-format off
PairStyle(ilp/water/2dm,PairILPWATER2DM);
// clang-format on
#else

#ifndef LMP_PAIR_ILP_WATER_2DM_H
#define LMP_PAIR_ILP_WATER_2DM_H

#include "pair_ilp_tmd.h"

namespace LAMMPS_NS {

class PairILPWATER2DM : virtual public PairILPTMD {
 public:
  PairILPWATER2DM(class LAMMPS *);

 protected:
  void settings(int, char **) override;

  /**************************************************************/
  /*       modulo operation with cycling around range           */

  inline int modulo(int k, int range)
  {
    if (k < 0) k += range;
    return k % range;
  }
  /**************************************************************/
};

}    // namespace LAMMPS_NS

#endif
#endif

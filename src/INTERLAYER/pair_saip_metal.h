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
PairStyle(saip/metal,PairSAIPMETAL);
// clang-format on
#else

#ifndef LMP_PAIR_SAIP_METAL_H
#define LMP_PAIR_SAIP_METAL_H

#include "pair_ilp_graphene_hbn.h"

namespace LAMMPS_NS {

class PairSAIPMETAL : virtual public PairILPGrapheneHBN {
 public:
  PairSAIPMETAL(class LAMMPS *);

 protected:
  void settings(int, char **) override;
  void calc_FRep(int, int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif

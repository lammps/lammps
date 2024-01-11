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
PairStyle(ilp/tmd,PairILPTMD);
// clang-format on
#else

#ifndef LMP_PAIR_ILP_TMD_H
#define LMP_PAIR_ILP_TMD_H

#include "pair_ilp_graphene_hbn.h"

namespace LAMMPS_NS {

class PairILPTMD : virtual public PairILPGrapheneHBN {
 public:
  PairILPTMD(class LAMMPS *);

 protected:
  void settings(int, char **) override;
  void ILP_neigh() override;
  void calc_normal() override;
  void calc_FRep(int, int) override;

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

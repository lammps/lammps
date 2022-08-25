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
PairStyle(srp/react,PairSRPREACT);
// clang-format on
#else

#ifndef LMP_PAIR_SRP_REACT_H
#define LMP_PAIR_SRP_REACT_H

#include "pair_srp.h"

namespace LAMMPS_NS {

class PairSRPREACT : public PairSRP {
 public:
  PairSRPREACT(class LAMMPS *);
  ~PairSRPREACT() override;
  void settings(int, char **) override;
  void init_style() override;

 private:
  char *idbreak;
  char *idcreate;
  bool bond_break, bond_create;
};
}    // namespace LAMMPS_NS
#endif
#endif

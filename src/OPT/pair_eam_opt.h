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
PairStyle(eam/opt,PairEAMOpt);
// clang-format on
#else

#ifndef LMP_PAIR_EAM_OPT_H
#define LMP_PAIR_EAM_OPT_H

#include "pair_eam.h"

namespace LAMMPS_NS {

// use virtual public since this class is parent in multiple inheritance

class PairEAMOpt : virtual public PairEAM {
 public:
  PairEAMOpt(class LAMMPS *);

  void compute(int, int) override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();
};

}    // namespace LAMMPS_NS

#endif
#endif

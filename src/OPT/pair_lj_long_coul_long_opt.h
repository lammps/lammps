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
PairStyle(lj/long/coul/long/opt,PairLJLongCoulLongOpt);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_LONG_COUL_LONG_OPT_H
#define LMP_PAIR_LJ_LONG_COUL_LONG_OPT_H

#include "pair_lj_long_coul_long.h"

namespace LAMMPS_NS {

class PairLJLongCoulLongOpt : public PairLJLongCoulLong {
 public:
  PairLJLongCoulLongOpt(class LAMMPS *);
  void compute(int, int) override;
  void compute_outer(int, int) override;

 protected:
  template <const int EVFLAG, const int EFLAG, const int NEWTON_PAIR, const int CTABLE,
            const int LJTABLE, const int ORDER1, const int ORDER6>
  void eval();

  template <const int EVFLAG, const int EFLAG, const int NEWTON_PAIR, const int CTABLE,
            const int LJTABLE, const int ORDER1, const int ORDER6>
  void eval_outer();
};

}    // namespace LAMMPS_NS

#endif
#endif

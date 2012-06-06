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

PairStyle(lj/cut/coul/long/opt,PairLJCutCoulLongOpt)

#else

#ifndef LMP_PAIR_LJ_CUT_COUL_LONG_OPT_H
#define LMP_PAIR_LJ_CUT_COUL_LONG_OPT_H

#include "pair_lj_cut_coul_long.h"

namespace LAMMPS_NS {

class PairLJCutCoulLongOpt : public PairLJCutCoulLong {
 public:
  PairLJCutCoulLongOpt(class LAMMPS *);
  virtual void compute(int, int);

 protected:
  template <const int EVFLAG, const int EFLAG,
            const int NEWTON_PAIR, const int CTABLE >
  void eval();
};

}

#endif
#endif

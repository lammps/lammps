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

PairStyle(lj/long/coul/long/opt,PairLJLongCoulLongOpt)

#else

#ifndef LMP_PAIR_LJ_LONG_COUL_LONG_OPT_H
#define LMP_PAIR_LJ_LONG_COUL_LONG_OPT_H

#include "pair_lj_long_coul_long.h"

namespace LAMMPS_NS {

class PairLJLongCoulLongOpt : public PairLJLongCoulLong {
 public:
  PairLJLongCoulLongOpt(class LAMMPS *);
  virtual void compute(int, int);
  virtual void compute_outer(int,int);

 protected:
  template <const int EVFLAG, const int EFLAG,
    const int NEWTON_PAIR, const int CTABLE, const int LJTABLE,
    const int ORDER1, const int ORDER6 >
  void eval();

  template <const int EVFLAG, const int EFLAG,
    const int NEWTON_PAIR, const int CTABLE, const int LJTABLE,
    const int ORDER1, const int ORDER6 >
  void eval_outer();
};

}

#endif
#endif

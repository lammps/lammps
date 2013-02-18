/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/cut/thr,PairLJCutThr)

#else

#ifndef LMP_PAIR_LJ_CUT_THR_H
#define LMP_PAIR_LJ_CUT_THR_H

#include "pair_lj_cut.h"

namespace LAMMPS_NS {

class PairLJCutThr : public PairLJCut {

 public:
  PairLJCutThr(class LAMMPS *);

  virtual void compute(int, int);
  virtual void init_style();

 private:
  template <int, int> void eval();

};
}

#endif
#endif

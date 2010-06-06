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

#ifdef PAIR_CLASS

PairStyle(yukawa/colloid/omp,PairYukawaColloidOMP)

#else

#ifndef LMP_PAIR_YUKAWA_COLLOID_OMP_H
#define LMP_PAIR_YUKAWA_COLLOID_OMP_H

#include "pair_yukawa_omp.h"

namespace LAMMPS_NS {

class PairYukawaColloidOMP : public PairYukawaOMP {
 public:
  PairYukawaColloidOMP(class LAMMPS *);
  ~PairYukawaColloidOMP() {}
  void compute(int, int);
  void init_style();
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double &);
 protected:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval_col();  
};

}

#endif
#endif

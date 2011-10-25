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

PairStyle(buck/coul/long/omp,PairBuckCoulLongOMP)

#else

#ifndef LMP_PAIR_BUCK_COUL_LONG_OMP_H
#define LMP_PAIR_BUCK_COUL_LONG_OMP_H

#include "pair_buck_coul_long.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairBuckCoulLongOMP : public PairBuckCoulLong, public ThrOMP {

 public:
  PairBuckCoulLongOMP(class LAMMPS *);

  virtual void compute(int, int);
  virtual double memory_usage();

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData * const thr);
};

}

#endif
#endif

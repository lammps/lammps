/* -*- c++ -*- -------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_OMP_H
#define LMP_PAIR_OMP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOMP : protected Pair {

 protected:
  double *eng_vdwl_thr;         // per thread accumulated vdw energy
  double *eng_coul_thr;         // per thread accumulated coulomb energies
  double **virial_thr;          // per thread virial
  double **eatom_thr;		// per thread per atom energy
  double ***vatom_thr;		// per thread per atom virial

  int maxeatom_thr, maxvatom_thr;
  
 public:
  PairOMP(class LAMMPS *);
  virtual ~PairOMP();
  virtual double memory_usage();

 protected:
  void ev_setup(int, int);
  void ev_reduce();
  void ev_tally_thr(int, int, int, int, double, double, double,
		    double, double, double, int);
};

}

#endif

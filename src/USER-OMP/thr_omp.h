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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifndef LMP_THR_OMP_H
#define LMP_THR_OMP_H

#include "pointers.h"

namespace LAMMPS_NS {

// forward declarations
class Pair;
class Dihedral;

class ThrOMP {

 protected:
  const int thr_style;
  enum {PAIR=1, BOND, ANGLE, DIHEDRAL, IMPROPER, KSPACE, FIX, COMPUTE};

  LAMMPS *lmp;           // reference to base lammps object.

  double *eng_vdwl_thr;  // per thread accumulated vdw energy
  double *eng_coul_thr;  // per thread accumulated coulomb energies
  double *eng_bond_thr;  // per thread accumlated bonded energy

  double **virial_thr;   // per thread virial
  double **eatom_thr;    // per thread per atom energy
  double ***vatom_thr;   // per thread per atom virial

  int maxeatom_thr, maxvatom_thr;
  
 public:
  ThrOMP(LAMMPS *, int);
  virtual ~ThrOMP();

  double memory_usage_thr();

 protected:
  // extra ev_tally work for threaded styles
  void ev_setup_thr(Pair *);
  void ev_setup_thr(Dihedral *);

  void ev_reduce_thr(Pair *);
  void ev_reduce_thr(Dihedral *);

 private:
  // internal method to be used by multiple ev_setup_thr() methods
  void ev_zero_acc_thr(int, int, int, int, int, int);

 protected:
  // threading adapted versions of the ev_tally infrastructure
  void ev_tally_thr(Pair *, int, int, int, int, double, double,
		    double, double, double, double, int);

 protected:
  // set loop range, thread id, and force array offset for threaded runs.
  double **loop_setup_thr(double **, int &, int &, int &, int, int, int);

  // reduce per thread forces into the first part of the force array
  void force_reduce_thr(double *, int, int, int);
};

}
#endif

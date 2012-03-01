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

PairStyle(lj/cut/coul/long/tip4p/omp,PairLJCutCoulLongTIP4POMP)

#else

#ifndef LMP_PAIR_LJ_CUT_COUL_LONG_TIP4P_OMP_H
#define LMP_PAIR_LJ_CUT_COUL_LONG_TIP4P_OMP_H

#include "pair_lj_cut_coul_long_tip4p.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairLJCutCoulLongTIP4POMP : public PairLJCutCoulLongTIP4P, public ThrOMP {

 public:
  PairLJCutCoulLongTIP4POMP(class LAMMPS *);
  virtual ~PairLJCutCoulLongTIP4POMP();

  virtual void compute(int, int);
  virtual double memory_usage();

 protected:
  // this is to cache m-shift corrected positions.
  int maxmpos;        // size of the following arrays
  int *h1idx, *h2idx; // local index of hydrogen atoms
  double **mpos;      // coordinates corrected for m-shift.
  void find_M_permissive(int, int &, int &, double *);

 private:
  template <int CFLAG, int EVFLAG, int EFLAG, int VFLAG>
  void eval(int ifrom, int ito, ThrData * const thr);
};

}

#endif
#endif

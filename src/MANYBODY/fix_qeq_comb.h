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

#ifdef FIX_CLASS

FixStyle(qeq/comb,FixQEQComb)

#else

#ifndef LMP_FIX_QEQ_COMB_H
#define LMP_FIX_QEQ_COMB_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixQEQComb : public Fix {
 public:
  FixQEQComb(class LAMMPS *, int, char **);
  ~FixQEQComb();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int,int,int);
  double memory_usage();

 private:
  int me,firstflag;
  double precision;
  int nlevels_respa;
  double ngroup;
  FILE *fp;

  class PairComb *comb;
  int nmax;
  double *qf,*q1,*q2;
};

}

#endif
#endif

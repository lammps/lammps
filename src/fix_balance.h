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

FixStyle(balance,FixBalance)

#else

#ifndef LMP_FIX_BALANCE_H
#define LMP_FIX_BALANCE_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixBalance : public Fix {
 public:
  FixBalance(class LAMMPS *, int, char **);
  ~FixBalance();
  int setmask();
  void init();
  void setup(int);
  void setup_pre_exchange();
  void pre_exchange();
  void pre_neighbor();
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  int nevery,lbstyle,nitermax;
  double thresh,stopthresh;
  char bstr[3];
  FILE *fp;

  double imbnow;                // current imbalance factor
  double imbprev;               // imbalance factor before last rebalancing
  double imbfinal;              // imbalance factor after last rebalancing
  int maxperproc;               // max atoms on any processor
  int itercount;                // iteration count of last call to Balance
  int kspace_flag;              // 1 if KSpace solver defined
  int pending;

  class Balance *balance;
  class Irregular *irregular;

  void rebalance();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix balance string is invalid

The string can only contain the characters "x", "y", or "z".

E: Fix balance string is invalid for 2d simulation

The string cannot contain the letter "z".

E: Cannot open fix balance output file

Self-explanatory.

*/

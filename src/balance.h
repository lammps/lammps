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

#ifdef COMMAND_CLASS

CommandStyle(balance,Balance)

#else

#ifndef LMP_BALANCE_H
#define LMP_BALANCE_H

#include "mpi.h"
#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Balance : protected Pointers {
 public:
  Balance(class LAMMPS *);
  ~Balance();
  void command(int, char **);
  void dynamic_setup(int, int, char *, double);
  int dynamic();
  double imbalance_nlocal(int &);

 private:
  int me,nprocs;
  int xflag,yflag,zflag,dflag;
  int nrepeat,niter;
  double thresh;
  char *bstr;
  double *user_xsplit,*user_ysplit,*user_zsplit;

  int *ops;
  int nops;
  double *splits[3];
  bigint *counts[3];
  double *cuts;
  bigint *onecount;
  MPI_Comm commslice[3];

  int *pcount,*allcount;

  FILE *fp;                      // for debug output
  bigint laststep;

  void dynamic_setup(char *);
  int dynamic_once();
  double imbalance_splits(int &);
  void stats(int, int, double *, bigint *);
  void adjust(int, bigint *, double *);
  int binary(double, int, double *);

  void dumpout(bigint);          // for debug output
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Balance command before simulation box is defined

The balance command cannot be used before a read_data, read_restart,
or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot balance in z dimension for 2d simulation

Self-explanatory.

E: Balance dynamic string is invalid

The string can only contain the characters "x", "y", or "z".

E: Balance dynamic string is invalid for 2d simulation

The string cannot contain the letter "z".

E: Lost atoms via balance: original %ld current %ld

This should not occur.  Report the problem to the developers.

E: Cannot open balance output file

This error message can only occur if debug options
are uncommented in src/balance.cpp.

*/

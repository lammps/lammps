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

  int count;
  int *pcount,*allcount;

  FILE *fp;                      // for debug output
  bigint laststep;

  double imbalance_splits(int &);
  double imbalance_nlocal(int &);
  void dynamic_setup(char *);
  void dynamic();
  void stats(int, int, double *, bigint *);
  void adjust(int, bigint *, double *);
  int binary(double, int, double *);

  void dumpout(bigint);          // for debug output
};

}

#endif
#endif

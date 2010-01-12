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

#ifdef COMPUTE_CLASS

ComputeStyle(reduce,ComputeReduce)

#else

#ifndef LMP_COMPUTE_REDUCE_H
#define LMP_COMPUTE_REDUCE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeReduce : public Compute {
 public:
  ComputeReduce(class LAMMPS *, int, char **);
  virtual ~ComputeReduce();
  void init();
  double compute_scalar();
  void compute_vector();
  double memory_usage();

 protected:
  int me;
  int mode,nvalues,iregion;
  int *which,*argindex,*flavor,*value2index;
  char **ids;
  double *onevec;
  int *replace,*indices,*owner;
  int index;

  int maxatom;
  double *varatom;

  struct Pair {
    double value;
    int proc;
  };
  Pair pairme,pairall;

  virtual double compute_one(int, int);
  virtual double count(int);
  void combine(double &, double, int);
};

}

#endif
#endif

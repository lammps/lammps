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

#ifndef COMPUTE_REDUCE_H
#define COMPUTE_REDUCE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeReduce : public Compute {
 public:
  ComputeReduce(class LAMMPS *, int, char **);
  ~ComputeReduce();
  void init();
  double compute_scalar();
  void compute_vector();
  double memory_usage();

 private:
  int mode,nvalues;
  int *which,*argindex,*value2index;
  char **ids;
  double *onevec;

  int maxatom;
  double *varatom;

  double compute_one(int);
  void combine(double &, double);
};

}

#endif

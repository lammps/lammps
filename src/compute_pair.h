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

ComputeStyle(pair,ComputePair)

#else

#ifndef LMP_COMPUTE_PAIR_H
#define LMP_COMPUTE_PAIR_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePair : public Compute {
 public:
  ComputePair(class LAMMPS *, int, char **);
  ~ComputePair();
  void init();
  double compute_scalar();
  void compute_vector();

 private:
  int evalue,npair;
  char *pstyle;
  class Pair *pair;
  double *one;
};

}

#endif
#endif

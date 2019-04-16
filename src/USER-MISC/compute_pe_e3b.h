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

#ifdef COMPUTE_CLASS

ComputeStyle(pe/e3b,ComputePEE3B)

#else

#ifndef LMP_COMPUTE_PEE3B_H
#define LMP_COMPUTE_PEE3B_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePEE3B : public Compute {

 public:
  ComputePEE3B(class LAMMPS *, int, char **);
  virtual ~ComputePEE3B();

  void init();

  double compute_scalar();
  void   compute_vector();

 private:
  class PairE3B *e3b;
};

}

#endif
#endif

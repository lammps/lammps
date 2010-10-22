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

ComputeStyle(ti,ComputeTI)

#else

#ifndef COMPUTE_TI_H
#define COMPUTE_TI_H

#include "compute.h"

namespace LAMMPS_NS {
  
class ComputeTI : public Compute {
 public: 
  ComputeTI(class LAMMPS *, int, char **);
  ~ComputeTI();
  void init();
  double compute_scalar();

 private:
  int nterms;
  int *which;
  int *ivar1,*ivar2;
  char **var1,**var2;
  class Pair **pptr;
  char **pstyle;
};

}

#endif
#endif

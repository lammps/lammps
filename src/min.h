/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef MIN_H
#define MIN_H

#include "lammps.h"

class Min : public LAMMPS {
 public:
  double einitial,efinal,eprevious;
  double gnorm2_init,gnorminf_init,gnorm2_final,gnorminf_final;
  int niter,neval;
  double dmin,dmax;
  int linestyle,lineiter;

  Min();
  virtual ~Min() {}
  virtual void init() = 0;
  virtual void run() = 0;
  virtual int memory_usage() {return 0;}

  void modify_params(int, char **);
};

#endif

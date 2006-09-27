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

#ifndef FIX_CENTRO_H
#define FIX_CENTRO_H

#include "fix.h"

class FixCentro : public Fix {
  friend class DumpCustom;

 public:
  FixCentro(int, char **);
  ~FixCentro();
  int setmask();
  void init();
  void dump();
  int memory_usage();

 private:
  int nmax,maxneigh;
  double *centro;
  double *distsq;
  int *nearest;

  void select(int, int, double *);
  void select2(int, int, double *, int *);
};

#endif

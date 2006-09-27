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

#ifndef RANMARS_H
#define RANMARS_H

#include "lammps.h"

class RanMars : public LAMMPS {
 public:
  RanMars(int);
  ~RanMars();
  double uniform();
  double gaussian();

 private:
  int seed,save;
  double second;
  double *u;
  int i97,j97;
  double c,cd,cm;
};

#endif

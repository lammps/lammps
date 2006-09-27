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

#ifndef SPECIAL_H
#define SPECIAL_H

#include "lammps.h"

class Special : public LAMMPS {
 public:
  Special();
  ~Special();
  void build();

 private:
  int me,nprocs;
  int **onetwo,**onethree,**onefour;

  void combine();
};

#endif

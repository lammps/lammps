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

#ifndef CREATE_ATOMS_H
#define CREATE_ATOMS_H

#include "lammps.h"

class CreateAtoms : public LAMMPS {
 public:
  CreateAtoms() {}
  ~CreateAtoms() {}
  void command(int, char **);

 private:
  int regionflag;
  double boxxhi,boxyhi,boxzhi;
  double subxlo,subxhi,subylo,subyhi,subzlo,subzhi;

  void add_atom(int, double, double, double);
};

#endif

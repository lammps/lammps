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


#ifndef WRITE_RESTART_H
#define WRITE_RESTART_H

#include "stdio.h"
#include "lammps.h"

class WriteRestart : public LAMMPS {
 public:
  WriteRestart();
  ~WriteRestart() {}
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  FILE *fp;
  double natoms;         // natoms (sum of nlocal) to write into file

  void header();
  void mass();
  void dipole();
  void force_fields();

  void write_int(int, int);
  void write_double(int, double);
  void write_char(int, char *);
};

#endif

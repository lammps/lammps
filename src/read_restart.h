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

#ifndef READ_RESTART_H
#define READ_RESTART_H

#include "stdio.h"
#include "lammps.h"

class ReadRestart : public LAMMPS {
 public:
  ReadRestart() {}
  ~ReadRestart() {}
  void command(int, char **);

 private:
  int me;
  FILE *fp;
  int nprocs_file;

  void file_search(char *, char *);
  void header();
  void mass();
  void dipole();
  void force_fields();

  int read_int();
  double read_double();
  char *read_char();
};

#endif

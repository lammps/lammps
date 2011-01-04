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

#ifdef COMMAND_CLASS

CommandStyle(write_restart,WriteRestart)

#else

#ifndef LMP_WRITE_RESTART_H
#define LMP_WRITE_RESTART_H

#include "stdio.h"
#include "pointers.h"
#include "lmptype.h"

namespace LAMMPS_NS {

class WriteRestart : protected Pointers {
 public:
  WriteRestart(class LAMMPS *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  FILE *fp;
  bigint natoms;         // natoms (sum of nlocal) to write into file

  void header();
  void type_arrays();
  void force_fields();

  void write_int(int, int);
  void write_double(int, double);
  void write_char(int, char *);
  void write_bigint(int, bigint);
};

}

#endif
#endif

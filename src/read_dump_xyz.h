/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributed by Timothy Sirk
------------------------------------------------------------------------- */

#ifndef LMP_READ_DUMP_XYZ_H
#define LMP_READ_DUMP_XYZ_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class ReadDumpXYZ : protected Pointers {
 public:
  ReadDumpXYZ(class LAMMPS *);
  ~ReadDumpXYZ();

  void file(FILE *);
  int read_time(bigint &);
  void skip();
  bigint read_header(double [3][3], int &, int, int, int *, char **,
                     int, int &, int &, int &, int &);
  void read_atoms(int, int, double **);

private:
  FILE *fp;                // pointer to file opened by caller
  char *line;              // line read from dump file
  bigint nstep;            // current (time) step number
  bigint natoms;           // current number of atoms
  bigint nid;              // current atom id.

  int *fieldindex;         // mapping of input fields to dump

  void read_lines(int);
};

}

#endif

/* ----------------------------------------------------------------------
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

#ifndef LMP_READ_DUMP_NATIVE_H
#define LMP_READ_DUMP_NATIVE_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class ReadDumpNative : protected Pointers {
 public:
  ReadDumpNative(class LAMMPS *);
  ~ReadDumpNative();

  void init(FILE *);
  void scan(bigint, int, int *, char **, int, bigint &, double [3][3], int &);
  void read(int, double **);

private:
  FILE *fp;                // pointer to file opened by caller
  char *line;              // line read from dump file
  int dimension;
  int triclinic;

  int nwords;              // # of per-atom columns in dump file
  char **words;            // ptrs to words in parsed per-atom line

  int nfield;              // # of fields to extract for each atom
  int *fieldindex;         // index into words for each field

  int find_label(const char *, int, char **);
  void read_lines(int);
};

}

#endif

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

#ifndef LMP_READER_H
#define LMP_READER_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Reader : protected Pointers {
 public:
  Reader(class LAMMPS *);
  virtual ~Reader() {}

  virtual int read_time(bigint &) = 0;
  virtual void skip() = 0;
  virtual bigint read_header(double [3][3], int &, int, int, int *, char **,
                             int, int &, int &, int &, int &) = 0;
  virtual void read_atoms(int, int, double **) = 0;

  void file(FILE *);

 protected:
  FILE *fp;                // pointer to file opened by caller
};

}

#endif

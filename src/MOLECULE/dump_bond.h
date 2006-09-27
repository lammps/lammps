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

#ifndef DUMP_BOND_H
#define DUMP_BOND_H

#include "dump.h"

class DumpBond : public Dump {
 public:
  DumpBond(int, char **);
  ~DumpBond() {}
  void init();

 private:
  int index;            // counter for bond output

  void write_header(int);
  int count();
  int pack();
  void write_data(int, double *);
};

#endif

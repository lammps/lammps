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

#ifndef FIX_COM_H
#define FIX_COM_H

#include "stdio.h"
#include "fix.h"

class FixCOM : public Fix {
 public:
  FixCOM(int, char **);
  ~FixCOM();
  int setmask();
  void init();
  void setup();
  void end_of_step();

 private:
  int me,first;
  FILE *fp;
  double masstotal;
};

#endif

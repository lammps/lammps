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

#ifndef FIX_FREEZE_H
#define FIX_FREEZE_H

#include "fix.h"

class FixFreeze : public Fix {
 public:
  FixFreeze(int, char **);
  ~FixFreeze() {}
  int setmask();
  void init();
  void setup();
  void post_force(int);
};

#endif

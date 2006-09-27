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

#ifndef FIX_NVE_GRAN_H
#define FIX_NVE_GRAN_H

#include "fix.h"

class FixNVEGran : public Fix {
 public:
  FixNVEGran(int, char **);
  ~FixNVEGran() {}
  int setmask();
  void init();
  void initial_integrate();
  void final_integrate();

 private:
  double dtv,dtf,dtfphi;
};

#endif

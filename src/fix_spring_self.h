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

#ifndef FIX_SPRING_SELF_H
#define FIX_SPRING_SELF_H

#include "fix.h"

class FixSpringSelf : public Fix {
 public:
  FixSpringSelf(int, char **);
  ~FixSpringSelf();
  int setmask();
  void init();
  void setup();
  void post_force(int);
  void post_force_respa(int, int, int);

  int memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 private:
  double k;
  double **xoriginal;         // original coords of atoms
  int nlevels_respa;
};

#endif

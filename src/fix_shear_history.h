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

#ifndef FIX_SHEAR_HISTORY_H
#define FIX_SHEAR_HISTORY_H

#include "fix.h"

class FixShearHistory : public Fix {
  friend class Neighbor;
  friend class FixPour;

 public:
  FixShearHistory(int, char **);
  ~FixShearHistory();
  int setmask();
  void init();
  void pre_exchange();

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
  int *npartner;                // # of touching partners of each atom
  int **partner;                // tags for the partners
  double ***shearpartner;       // 3 shear values with the partner
};

#endif

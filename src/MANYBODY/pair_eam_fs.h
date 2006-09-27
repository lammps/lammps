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

#ifndef PAIR_EAM_FS_H
#define PAIR_EAM_FS_H

#include "pair_eam.h"

class PairEAMFS : public PairEAM {
 public:
  PairEAMFS();
  ~PairEAMFS();
  void compute(int, int);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void single(int, int, int, int, double, double, double, int, One &);

 private:
  double ***rhor_fs;
  double ***rhor_fs_0,***rhor_fs_1,***rhor_fs_2,***rhor_fs_3;
  double ***rhor_fs_4,***rhor_fs_5,***rhor_fs_6;

  int read_setfl(char *, int, int);
  void store_setfl();
  void interpolate();
  void interpolate_deallocate();
};

#endif

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_BONDSWAP_H
#define FIX_BONDSWAP_H

#include "fix.h"
#include "pair.h"

namespace LAMMPS_NS {

class FixBondSwap : public Fix {
 public:
  FixBondSwap(class LAMMPS *, int, char **);
  ~FixBondSwap();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void pre_neighbor();
  int modify_param(int, char **);
  double compute_vector(int);
  double memory_usage();

 private:
  double fraction,cutsq;
  int nmax,tflag;
  int *alist;
  int naccept,foursome;
  int angleflag;
  char *id_temp;
  int *type;
  double **x;
  
  class NeighList *list;
  class Compute *temperature;
  class RanMars *random;

  double dist_rsq(int, int);
  double pair_eng(int, int);
  double bond_eng(int, int, int);
  double angle_eng(int, int, int, int);
};

}

#endif

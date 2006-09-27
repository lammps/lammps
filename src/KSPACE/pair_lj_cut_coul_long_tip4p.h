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

#ifndef PAIR_LJ_CUT_COUL_LONG_TIP4P_H
#define PAIR_LJ_CUT_COUL_LONG_TIP4P_H

#include "pair_lj_cut_coul_long.h"

class PairLJCutCoulLongTIP4P : public PairLJCutCoulLong {
  friend class PPPM; 
  
 public:
  PairLJCutCoulLongTIP4P();
  void compute(int, int);
  void settings(int, char **);
  void init_style();
  void write_restart_settings(FILE *fp);
  void read_restart_settings(FILE *fp);

 private:
  int typeH,typeO;             // atom types of TIP4P water H and O atoms
  int typeA,typeB;             // angle and bond types of TIP4P water
  double qdist;                // distance from O site to negative charge
  double alpha;                // geometric constraint parameter for TIP4P
  double **tf;
  double tvirial[6];

  void find_M(int, int &, int &, double *);
};

#endif

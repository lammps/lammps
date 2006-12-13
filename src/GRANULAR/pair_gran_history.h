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

#ifndef PAIR_GRAN_HISTORY_H
#define PAIR_GRAN_HISTORY_H

#include "pair.h"

class PairGranHistory : public Pair {
  friend class Neighbor;
  friend class FixWallGran;
  friend class FixGranDiag;

 public:
  PairGranHistory();
  ~PairGranHistory();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 protected:
  double xkk,xkkt,xmu;
  int dampflag;
  double gamman;
  double dt,gamman_dl,gammas_dl;
  int ifix_history,freeze_group_bit;

  void allocate();
};

#endif

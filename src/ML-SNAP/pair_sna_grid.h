/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(sna/grid, PairSNAGrid);
// clang-format on
#else

#ifndef LMP_PAIR_SNA_GRID_H
#define LMP_PAIR_SNA_GRID_H

#include "pair_grid.h"

namespace LAMMPS_NS {

class PairSNAGrid : public PairGrid {
 public:
  PairSNAGrid(class LAMMPS *);
  ~PairSNAGrid();

  void init_style();
  void init_list(int, class NeighList *);
  void settings(int, char **);
  void compute(int, int);
  double memory_usage();

 private:
  int ncoeff;
  class NeighList *list;
  double rcutfac;
  double *radelem;
  double *wjelem;
  class SNA *snaptr;
  int quadraticflag;
  int twojmax, switchflag, bzeroflag, bnormflag;
  int chemflag, wselfallflag;
  int switchinnerflag;
  double *sinnerelem;
  double *dinnerelem;
  double rfac0, rmin0;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute sna/grid/local requires a pair style be defined

Self-explanatory.

E: Compute sna/grid/local cutoff is longer than pairwise cutoff

Self-explanatory.

W: More than one compute sna/grid/local

Self-explanatory.

*/

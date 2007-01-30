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

#ifndef PAIR_HYBRID_H
#define PAIR_HYBRID_H

#include "stdio.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairHybrid : public Pair {
 public:
  int nstyles;                  // # of different pair styles
  Pair **styles;                // class list for each Pair style
  char **keywords;              // sub-style name for each Pair style

  PairHybrid(class LAMMPS *);
  ~PairHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void single(int, int, int, int, double, double, double, int, One &);
  void single_embed(int, int, double &);
  void modify_params(int narg, char **arg);
  int memory_usage();

 private:
  int **map;                    // which style each itype,jtype points to

  int *nnlist;                  // # of half neighs in sub-style neigh lists
  int *maxneigh;                // max # of neighs sub-style lists can store
  int **nlist;                  // half neigh list for each sub-style

  int *nnlist_full;             // # of full neighs in sub-style neigh lists
  int *maxneigh_full;           // max # of neighs sub-style lists can store
  int **nlist_full;             // full neigh list for each sub-style

  int ***firstneigh;            // each sub-style's per-atom firstneigh
  int **numneigh;               // each sub-style's per-atom numneigh
  int ***firstneigh_full;       // each sub-style's per-atom firstneigh_full
  int **numneigh_full;          // each sub-style's per-atom numneigh_full
  int maxlocal;                 // max length of each ss's firstneigh,numneigh

  void allocate();
};

}

#endif

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
  int nstyles;                  // # of different sub-styles
  Pair **styles;                // list of Pair style classes
  char **keywords;              // style name of each Pair style

  PairHybrid(class LAMMPS *);
  ~PairHybrid();
  void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void modify_params(int narg, char **arg);
  double memory_usage();

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);
  void *extract(char *);
  void reset_dt();

 protected:
  int **nmap;                   // # of sub-styles itype,jtype points to
  int ***map;                   // list of sub-styles itype,jtype points to

  void allocate();
  virtual void modify_requests();
};

}

#endif

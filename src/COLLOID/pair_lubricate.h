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

#ifndef PAIR_LUBRICATE_H
#define PAIR_LUBRICATE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLubricate : public Pair {
 public:
  PairLubricate(class LAMMPS *);
  ~PairLubricate();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 protected:
  double cut_inner_global,cut_global;
  double t_target,mu;
  int flag1,flag2,flag3,flag4;
  int seed,omega_flag;
  double **cut_inner,**cut;

  class RanMars *random;

  void allocate();
};

}

#endif

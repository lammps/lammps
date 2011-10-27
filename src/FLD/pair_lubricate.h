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

#ifdef PAIR_CLASS

PairStyle(lubricate,PairLubricate)

#else

#ifndef LMP_PAIR_LUBRICATE_H
#define LMP_PAIR_LUBRICATE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLubricate : public Pair {
 public:
  PairLubricate(class LAMMPS *);
  virtual ~PairLubricate();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  virtual void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  int pre_adapt(char *, int, int, int, int);
  void adapt(int, int, int, int, int, double);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 protected:
  double mu,cut_inner_global,cut_global;
  int flaglog,flagfld,shearing;
  double Ef[3][3];
  double R0,RT0,RS0;
  double **cut_inner,**cut;

  void allocate();
};

}

#endif
#endif

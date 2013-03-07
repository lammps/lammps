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

PairStyle(lj/cut/tip4p/cut,PairLJCutTIP4PCut)

#else

#ifndef LMP_PAIR_LJ_CUT_TIP4P_CUT_H
#define LMP_PAIR_LJ_CUT_TIP4P_CUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutTIP4PCut : public Pair {
 public:
  PairLJCutTIP4PCut(class LAMMPS *);
  virtual ~PairLJCutTIP4PCut();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double memory_usage();

 protected:
  double cut_lj_global,cut_coul_global;
  double cut_coul,cut_coulsq;
  double cut_coulsqplus;       // extended value for cut_coulsq
  double **cut_lj,**cut_ljsq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;

  int typeH,typeO;             // atom types of TIP4P water H and O atoms
  int typeA,typeB;             // angle and bond types of TIP4P water
  double alpha;                // geometric constraint parameter for TIP4P
  double qdist;

  int nmax;                    // info on off-oxygen charge sites
  int **hneigh;                // 0,1 = indices of 2 H associated with O
                               // 2 = 0 if site loc not yet computed, 1 if yes
  double **newsite;            // locations of charge sites

  void allocate();
  void compute_newsite(double *, double *, double *, double *);
};
}

#endif
#endif

/* -*- c++ -*- ----------------------------------------------------------
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

PairStyle(oxdna/stk,PairOxdnaStk)
PairStyle(oxdna2/stk,PairOxdnaStk)

#else

#ifndef LMP_PAIR_OXDNA_STK_H
#define LMP_PAIR_OXDNA_STK_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxdnaStk : public Pair {
 public:
  PairOxdnaStk(class LAMMPS *);
  virtual ~PairOxdnaStk();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  void *extract(const char *, int &);

 protected:
  // stacking interaction
  double stacking_strength(double, double, double);
  double **epsilon_st, **a_st, **cut_st_0, **cut_st_c;
  double **cut_st_lo, **cut_st_hi;
  double **cut_st_lc, **cut_st_hc, **b_st_lo, **b_st_hi, **shift_st;
  double **cutsq_st_hc;
  double **a_st4, **theta_st4_0, **dtheta_st4_ast;
  double **b_st4, **dtheta_st4_c;
  double **a_st5, **theta_st5_0, **dtheta_st5_ast;
  double **b_st5, **dtheta_st5_c;
  double **a_st6, **theta_st6_0, **dtheta_st6_ast;
  double **b_st6, **dtheta_st6_c;
  double **a_st1, **cosphi_st1_ast, **b_st1, **cosphi_st1_c;
  double **a_st2, **cosphi_st2_ast, **b_st2, **cosphi_st2_c;

  int seqdepflag;

  virtual void allocate();
  void ev_tally_xyz(int, int, int, int, double, double, double, double, double, double, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/

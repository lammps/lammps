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

PairStyle(oxrna2/stk,PairOxrna2Stk)

#else

#ifndef LMP_PAIR_OXRNA2_STK_H
#define LMP_PAIR_OXRNA2_STK_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxrna2Stk : public Pair {
 public:
  PairOxrna2Stk(class LAMMPS *);
  virtual ~PairOxrna2Stk();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
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
  double eta_st[4][4];
  double stacking_strength(double, double, double);
  double **epsilon_st, **a_st, **cut_st_0, **cut_st_c;
  double **cut_st_lo, **cut_st_hi;
  double **cut_st_lc, **cut_st_hc, **b_st_lo, **b_st_hi, **shift_st;
  double **cutsq_st_hc;
  double **a_st5, **theta_st5_0, **dtheta_st5_ast;
  double **b_st5, **dtheta_st5_c;
  double **a_st6, **theta_st6_0, **dtheta_st6_ast;
  double **b_st6, **dtheta_st6_c;
  double **a_st9, **theta_st9_0, **dtheta_st9_ast;
  double **b_st9, **dtheta_st9_c;
  double **a_st10, **theta_st10_0, **dtheta_st10_ast;
  double **b_st10, **dtheta_st10_c;
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

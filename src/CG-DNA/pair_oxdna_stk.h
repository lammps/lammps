/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(oxdna/stk,PairOxdnaStk);
PairStyle(oxdna2/stk,PairOxdnaStk);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA_STK_H
#define LMP_PAIR_OXDNA_STK_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxdnaStk : public Pair {
 public:
  PairOxdnaStk(class LAMMPS *);
  ~PairOxdnaStk() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void init_list(int, class NeighList *) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  void *extract(const char *, int &) override;

 protected:
  // stacking interaction
  double eta_st[4][4];
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
  double **nx_xtrct, **ny_xtrct, **nz_xtrct;    // per-atom arrays for local unit vectors
  int seqdepflag;

  virtual void allocate();
  void ev_tally_xyz(int, int, int, int, double, double, double, double, double, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif

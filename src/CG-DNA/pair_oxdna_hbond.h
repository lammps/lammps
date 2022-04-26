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
PairStyle(oxdna/hbond,PairOxdnaHbond);
PairStyle(oxdna2/hbond,PairOxdnaHbond);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA_HBOND_H
#define LMP_PAIR_OXDNA_HBOND_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxdnaHbond : public Pair {
 public:
  PairOxdnaHbond(class LAMMPS *);
  ~PairOxdnaHbond() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
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
  // h-bonding interaction
  double alpha_hb[4][4];
  double **epsilon_hb, **a_hb, **cut_hb_0, **cut_hb_c, **cut_hb_lo, **cut_hb_hi;
  double **cut_hb_lc, **cut_hb_hc, **b_hb_lo, **b_hb_hi, **shift_hb;
  double **cutsq_hb_hc;
  double **a_hb1, **theta_hb1_0, **dtheta_hb1_ast;
  double **b_hb1, **dtheta_hb1_c;
  double **a_hb2, **theta_hb2_0, **dtheta_hb2_ast;
  double **b_hb2, **dtheta_hb2_c;
  double **a_hb3, **theta_hb3_0, **dtheta_hb3_ast;
  double **b_hb3, **dtheta_hb3_c;
  double **a_hb4, **theta_hb4_0, **dtheta_hb4_ast;
  double **b_hb4, **dtheta_hb4_c;
  double **a_hb7, **theta_hb7_0, **dtheta_hb7_ast;
  double **b_hb7, **dtheta_hb7_c;
  double **a_hb8, **theta_hb8_0, **dtheta_hb8_ast;
  double **b_hb8, **dtheta_hb8_c;
  double **nx_xtrct, **ny_xtrct, **nz_xtrct;    // per-atom arrays for local unit vectors
  int seqdepflag;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

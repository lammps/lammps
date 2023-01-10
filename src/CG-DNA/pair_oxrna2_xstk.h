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
PairStyle(oxrna2/xstk,PairOxrna2Xstk);
// clang-format on
#else

#ifndef LMP_PAIR_OXRNA2_XSTK_H
#define LMP_PAIR_OXRNA2_XSTK_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxrna2Xstk : public Pair {
 public:
  PairOxrna2Xstk(class LAMMPS *);
  ~PairOxrna2Xstk() override;
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
  // cross-stacking interaction
  double **k_xst, **cut_xst_0, **cut_xst_c, **cut_xst_lo, **cut_xst_hi;
  double **cut_xst_lc, **cut_xst_hc, **b_xst_lo, **b_xst_hi;
  double **cutsq_xst_hc;
  double **a_xst1, **theta_xst1_0, **dtheta_xst1_ast;
  double **b_xst1, **dtheta_xst1_c;
  double **a_xst2, **theta_xst2_0, **dtheta_xst2_ast;
  double **b_xst2, **dtheta_xst2_c;
  double **a_xst3, **theta_xst3_0, **dtheta_xst3_ast;
  double **b_xst3, **dtheta_xst3_c;
  double **a_xst7, **theta_xst7_0, **dtheta_xst7_ast;
  double **b_xst7, **dtheta_xst7_c;
  double **a_xst8, **theta_xst8_0, **dtheta_xst8_ast;
  double **b_xst8, **dtheta_xst8_c;
  double **nx_xtrct, **nz_xtrct;    // per-atom arrays for local unit vectors

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

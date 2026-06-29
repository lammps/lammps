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
PairStyle(oxdna3/xstk,PairOxdna3Xstk);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA3_XSTK_H
#define LMP_PAIR_OXDNA3_XSTK_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxdna3Xstk : public Pair {
 public:
  PairOxdna3Xstk(class LAMMPS *);
  ~PairOxdna3Xstk() override;
  void compute_base_site(int, double *, double *, double *, double *) const;
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
  void *extract(const char *, int &) override;

 protected:
  // cross-stacking interaction
  double **k_xst;
  double ****cut_xst_0_33, ****cut_xst_c_33, ****cut_xst_lo_33, ****cut_xst_hi_33;
  double ****cut_xst_lc_33, ****cut_xst_hc_33, ****cutsq_xst_hc_33;
  double ****cut_xst_0_55, ****cut_xst_c_55, ****cut_xst_lo_55, ****cut_xst_hi_55;
  double ****cut_xst_lc_55, ****cut_xst_hc_55, ****cutsq_xst_hc_55;
  double **b_xst_lo, **b_xst_hi;
  double **a_xst1, **theta_xst1_0, **dtheta_xst1_ast;
  double **b_xst1, **dtheta_xst1_c;
  double **a_xst2, **theta_xst2_0, **dtheta_xst2_ast;
  double **b_xst2, **dtheta_xst2_c;
  double **a_xst3, **theta_xst3_0, **dtheta_xst3_ast;
  double **b_xst3, **dtheta_xst3_c;
  double ****a_xst4_33, ****theta_xst4_0_33, ****dtheta_xst4_ast_33;
  double ****b_xst4_33, ****dtheta_xst4_c_33;
  double ****a_xst4_55, ****theta_xst4_0_55, ****dtheta_xst4_ast_55;
  double ****b_xst4_55, ****dtheta_xst4_c_55;
  double **a_xst7, **theta_xst7_0_33, **theta_xst7_0_55, **dtheta_xst7_ast;
  double **b_xst7, **dtheta_xst7_c;
  double **a_xst8, **theta_xst8_0_33, **theta_xst8_0_55, **dtheta_xst8_ast;
  double **b_xst8, **dtheta_xst8_c;
  double **nxyz_xtrct;    // per-atom arrays for local unit vectors

  virtual void allocate();

  class FixOxdnaLRF *fix_lrf;    // ptr to oxdna/lrf fix
};

}    // namespace LAMMPS_NS

#endif
#endif

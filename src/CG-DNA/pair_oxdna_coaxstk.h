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
PairStyle(oxdna/coaxstk,PairOxdnaCoaxstk);
PairStyle(oxrna2/coaxstk,PairOxdnaCoaxstk);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA_COAXSTK_H
#define LMP_PAIR_OXDNA_COAXSTK_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxdnaCoaxstk : public Pair {
 public:
  PairOxdnaCoaxstk(class LAMMPS *);
  ~PairOxdnaCoaxstk() override;
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
  // coaxial stacking interaction
  double **k_cxst, **cut_cxst_0, **cut_cxst_c, **cut_cxst_lo, **cut_cxst_hi;
  double **cut_cxst_lc, **cut_cxst_hc, **b_cxst_lo, **b_cxst_hi;
  double **cutsq_cxst_hc;
  double **a_cxst1, **theta_cxst1_0, **dtheta_cxst1_ast;
  double **b_cxst1, **dtheta_cxst1_c;
  double **a_cxst4, **theta_cxst4_0, **dtheta_cxst4_ast;
  double **b_cxst4, **dtheta_cxst4_c;
  double **a_cxst5, **theta_cxst5_0, **dtheta_cxst5_ast;
  double **b_cxst5, **dtheta_cxst5_c;
  double **a_cxst6, **theta_cxst6_0, **dtheta_cxst6_ast;
  double **b_cxst6, **dtheta_cxst6_c;
  double **a_cxst3p, **cosphi_cxst3p_ast, **b_cxst3p, **cosphi_cxst3p_c;
  double **a_cxst4p, **cosphi_cxst4p_ast, **b_cxst4p, **cosphi_cxst4p_c;
  double **nx_xtrct, **ny_xtrct, **nz_xtrct;    // per-atom arrays for local unit vectors

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

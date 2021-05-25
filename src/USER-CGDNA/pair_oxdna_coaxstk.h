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
  virtual ~PairOxdnaCoaxstk();
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

  virtual void allocate();
};

}    // namespace LAMMPS_NS

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

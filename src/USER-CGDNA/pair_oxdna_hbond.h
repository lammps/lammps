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

PairStyle(oxdna/hbond,PairOxdnaHbond)
PairStyle(oxdna2/hbond,PairOxdnaHbond)

#else

#ifndef LMP_PAIR_OXDNA_HBOND_H
#define LMP_PAIR_OXDNA_HBOND_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxdnaHbond : public Pair {
 public:
  PairOxdnaHbond(class LAMMPS *);
  virtual ~PairOxdnaHbond();
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
  // h-bonding interaction
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

  int seqdepflag;

  virtual void allocate();
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

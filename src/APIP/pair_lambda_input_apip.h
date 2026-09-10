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
/* ----------------------------------------------------------------------
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lambda/input/apip,PairLambdaInputAPIP);
// clang-format on
#else

#ifndef LMP_PAIR_LAMBDA_INPUT_APIP_H
#define LMP_PAIR_LAMBDA_INPUT_APIP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLambdaInputAPIP : public Pair {
  friend class FixLambdaAPIP;

 public:
  PairLambdaInputAPIP(class LAMMPS *);
  ~PairLambdaInputAPIP() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;

 protected:
  class FixLambdaAPIP *fix_lambda;    // ptr to fix lambda to store the calculated lambda_input

  // pro forma pair style variables
  double cut_global;
  double **cut;

  double timer;            // accumulated compute time
  double time_per_atom;    // compute time of one calculaiton
  int n_calculations;      // number of accumulated calculations
  int ignore_group_bit;    // groupbit of the group that must not be calculated

  void allocate();
  void calculate_time_per_atom();
  virtual int calculate_lambda_input();
};

}    // namespace LAMMPS_NS

#endif
#endif

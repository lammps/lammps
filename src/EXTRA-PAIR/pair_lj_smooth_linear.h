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
PairStyle(lj/smooth/linear,PairLJSmoothLinear);
PairStyle(lj/sf,PairLJSmoothLinear);
// clang-format on
#else

#ifndef PAIR_LJ_SMOOTH_LINEAR_H
#define PAIR_LJ_SMOOTH_LINEAR_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJSmoothLinear : public Pair {
 public:
  PairLJSmoothLinear(class LAMMPS *);
  ~PairLJSmoothLinear() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  double single_hessian(int, int, int, int, double, double[3], double, double, double &,
                        double[6]) override;

 protected:
  double cut_global;
  double **cut;
  double **epsilon, **sigma;
  double **ljcut, **dljcut;
  double **lj1, **lj2, **lj3, **lj4;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

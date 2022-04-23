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
PairStyle(hbond/dreiding/lj,PairHbondDreidingLJ);
// clang-format on
#else

#ifndef LMP_PAIR_HBOND_DREIDING_LJ_H
#define LMP_PAIR_HBOND_DREIDING_LJ_H

#include "pair.h"

namespace LAMMPS_NS {

class PairHbondDreidingLJ : public Pair {
 public:
  PairHbondDreidingLJ(class LAMMPS *);
  ~PairHbondDreidingLJ() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  double cut_inner_global, cut_outer_global, cut_angle_global;
  int ap_global;

  struct Param {
    double epsilon, sigma;
    double lj1, lj2, lj3, lj4;
    double d0, alpha, r0;
    double morse1;
    double denom_vdw;
    double cut_inner, cut_outer, cut_innersq, cut_outersq, cut_angle, offset;
    int ap;
  };

  Param *params;    // parameter set for an I-J-K interaction
  int nparams;      // number of parameters read
  int maxparam;

  int *donor;           // 1 if this type is ever a donor, else 0
  int *acceptor;        // 1 if this type is ever an acceptor, else 0
  int ***type2param;    // mapping from D,A,H to params, -1 if no map

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

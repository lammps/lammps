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
PairStyle(lj/cut/coul/cut/soft,PairLJCutCoulCutSoft);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_CUT_SOFT_H
#define LMP_PAIR_LJ_CUT_COUL_CUT_SOFT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutCoulCutSoft : public Pair {
 public:
  PairLJCutCoulCutSoft(class LAMMPS *);
  ~PairLJCutCoulCutSoft() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

 protected:
  double cut_lj_global, cut_coul_global;
  double **cut_lj, **cut_ljsq;
  double **cut_coul, **cut_coulsq;
  double **epsilon, **sigma, **lambda;
  double nlambda, alphalj, alphac;
  double **lj1, **lj2, **lj3, **lj4, **offset;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

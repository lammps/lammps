/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/dipole/long,PairLJCutDipoleLong);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_DIPOLE_LONG_H
#define LMP_PAIR_LJ_CUT_DIPOLE_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutDipoleLong : public Pair {
 public:
  double cut_coul;
  double **sigma;

  PairLJCutDipoleLong(class LAMMPS *);
  ~PairLJCutDipoleLong() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void *extract(const char *, int &) override;

 protected:
  double cut_lj_global;
  double **cut_lj, **cut_ljsq;
  double cut_coulsq;
  double **epsilon;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  double g_ewald;
  int ewald_order;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

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
PairStyle(lj/long/dipole/long,PairLJLongDipoleLong);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_LONG_DIPOLE_LONG_H
#define LMP_PAIR_LJ_LONG_DIPOLE_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJLongDipoleLong : public Pair {
 public:
  double cut_coul;

  PairLJLongDipoleLong(class LAMMPS *);
  ~PairLJLongDipoleLong() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;

  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void *extract(const char *, int &) override;

 protected:
  double cut_lj_global;
  double **cut_lj, **cut_lj_read, **cut_ljsq;
  double cut_coulsq;
  double **epsilon_read, **epsilon, **sigma_read, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  double *cut_respa;
  double g_ewald;
  int ewald_order, ewald_off;

  void options(char **arg, int order);
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

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
PairStyle(buck/long/coul/long,PairBuckLongCoulLong);
// clang-format on
#else

#ifndef LMP_PAIR_BUCK_LONG_COUL_LONG_H
#define LMP_PAIR_BUCK_LONG_COUL_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBuckLongCoulLong : public Pair {
 public:
  double cut_coul;

  PairBuckLongCoulLong(class LAMMPS *);
  ~PairBuckLongCoulLong() override;
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

  void compute_inner() override;
  void compute_middle() override;
  void compute_outer(int, int) override;

 protected:
  double cut_buck_global;
  double **cut_buck, **cut_buck_read, **cut_bucksq;
  double cut_coulsq;
  double **buck_a_read, **buck_a, **buck_c_read, **buck_c;
  double **buck1, **buck2, **buck_rho_read, **buck_rho, **rhoinv, **offset;
  double *cut_respa;
  double g_ewald;
  double g_ewald_6;
  int ewald_order, ewald_off;

  void options(char **arg, int order);
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

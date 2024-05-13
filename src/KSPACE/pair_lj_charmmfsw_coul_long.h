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
PairStyle(lj/charmmfsw/coul/long,PairLJCharmmfswCoulLong);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CHARMMFSW_COUL_LONG_H
#define LMP_PAIR_LJ_CHARMMFSW_COUL_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCharmmfswCoulLong : public Pair {
 public:
  PairLJCharmmfswCoulLong(class LAMMPS *);
  ~PairLJCharmmfswCoulLong() override;

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

  void compute_inner() override;
  void compute_middle() override;
  void compute_outer(int, int) override;
  void *extract(const char *, int &) override;

 protected:
  int implicit;
  int dihedflag;

  double cut_lj_inner, cut_lj, cut_ljinv, cut_lj_innerinv;
  double cut_lj_innersq, cut_ljsq;
  double cut_lj3inv, cut_lj_inner3inv, cut_lj3, cut_lj_inner3;
  double cut_lj6inv, cut_lj_inner6inv, cut_lj6, cut_lj_inner6;
  double cut_coul, cut_coulsq;
  double cut_bothsq;
  double denom_lj, denom_lj12, denom_lj6;
  double **epsilon, **sigma, **eps14, **sigma14;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  double **lj14_1, **lj14_2, **lj14_3, **lj14_4;
  double *cut_respa;
  double g_ewald;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

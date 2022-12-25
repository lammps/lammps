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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/spica/coul/long,PairLJSPICACoulLong);
PairStyle(lj/sdk/coul/long,PairLJSPICACoulLong);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_SPICA_COUL_LONG_H
#define LMP_PAIR_LJ_SPICA_COUL_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJSPICACoulLong : public Pair {
 public:
  PairLJSPICACoulLong(class LAMMPS *);
  ~PairLJSPICACoulLong() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;
  double memory_usage() override;

 protected:
  double **cut_lj, **cut_ljsq;
  double cut_coul, cut_coulsq;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  int **lj_type;

  // cutoff and offset for minimum of LJ potential
  // to be used in SPICA angle potential, which
  // uses only the repulsive part of the potential

  double **rminsq, **emin;

  double cut_lj_global;
  double g_ewald;

  void allocate();

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();
};

}    // namespace LAMMPS_NS

#endif
#endif

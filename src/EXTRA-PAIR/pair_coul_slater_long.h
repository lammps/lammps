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
PairStyle(coul/slater/long,PairCoulSlaterLong);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_SLATER_LONG_H
#define LMP_PAIR_COUL_SLATER_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulSlaterLong : public Pair {
 public:
  PairCoulSlaterLong(class LAMMPS *);
  ~PairCoulSlaterLong() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

 protected:
  double cut_coul, cut_coulsq, lamda, g_ewald;
  double **scale;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

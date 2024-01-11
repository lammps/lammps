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
PairStyle(colloid,PairColloid);
// clang-format on
#else

#ifndef LMP_PAIR_COLLOID_H
#define LMP_PAIR_COLLOID_H

#include "pair.h"

namespace LAMMPS_NS {

class PairColloid : public Pair {
 public:
  PairColloid(class LAMMPS *);
  ~PairColloid() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  enum { SMALL_SMALL, SMALL_LARGE, LARGE_LARGE };

  double cut_global;
  double **cut;
  double **a12, **d1, **d2, **diameter, **a1, **a2, **offset;
  double **sigma, **sigma3, **sigma6;
  double **lj1, **lj2, **lj3, **lj4;
  int **form;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

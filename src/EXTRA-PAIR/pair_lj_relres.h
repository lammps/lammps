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
PairStyle(lj/relres,PairLJRelRes);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_RELRES_H
#define LMP_PAIR_LJ_RELRES_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJRelRes : public Pair {
 public:
  PairLJRelRes(class LAMMPS *);
  ~PairLJRelRes() override;
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
  double cut_inner_global, cut_global, cutf_inner_global, cutf_global;
  double **cut, **cut_inner, **cut_inner_sq, **cutf, **cutfsq, **cutf_inner, **cutf_inner_sq;
  double **epsilon, **sigma;
  double **epsilonf, **sigmaf;
  double **lj1, **lj2, **lj3, **lj4;
  double **ljf1, **ljf2, **ljf3, **ljf4;
  double **ljsw0, **ljsw1, **ljsw2, **ljsw3, **ljsw4;
  double **ljswf0, **ljswf1, **ljswf2, **ljswf3, **ljswf4;
  double **offset, **offsetsp, **offsetsm;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

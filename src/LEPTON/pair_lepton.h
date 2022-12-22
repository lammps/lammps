/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Pair zero is a dummy pair interaction useful for requiring a
   force cutoff distance in the absence of pair-interactions or
   with hybrid/overlay if a larger force cutoff distance is required.

   This can be used in conjunction with bond/create to create bonds
   that are longer than the cutoff of a given force field, or to
   calculate radial distribution functions for models without
   pair interactions.

------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lepton,PairLepton);
// clang-format on
#else

#ifndef LMP_PAIR_LEPTON_H
#define LMP_PAIR_LEPTON_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLepton : public Pair {
 public:
  PairLepton(class LAMMPS *);
  ~PairLepton() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  std::vector<std::string> expressions;
  double **cut;
  int **type2expression;
  double cut_global;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

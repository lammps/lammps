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
PairStyle(hybrid/scaled,PairHybridScaled);
PairStyle(hybrid/scaled/omp,PairHybridScaled);
// clang-format on
#else

#ifndef LMP_PAIR_HYBRID_SCALED_H
#define LMP_PAIR_HYBRID_SCALED_H

#include "pair_hybrid.h"

#include <string>
#include <vector>

namespace LAMMPS_NS {

class PairHybridScaled : public PairHybrid {
 public:
  PairHybridScaled(class LAMMPS *);
  ~PairHybridScaled() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void born_matrix(int, int, int, int, double, double, double, double &, double &) override;

  void init_svector() override;
  void copy_svector(int, int) override;

 protected:
  double **fsum, **tsum;
  double *scaleval;
  int *scaleidx;
  std::vector<std::string> scalevars;
  int nmaxfsum;
};

}    // namespace LAMMPS_NS

#endif
#endif

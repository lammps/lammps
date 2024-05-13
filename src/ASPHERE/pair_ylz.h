/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(ylz,PairYLZ);
// clang-format on

#else

#ifndef LMP_PAIR_YLZ_H
#define LMP_PAIR_YLZ_H

#include "pair.h"

namespace LAMMPS_NS {

class PairYLZ : public Pair {
 public:
  PairYLZ(LAMMPS *lmp);
  ~PairYLZ() override;
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

 protected:
  double cut_global;
  double **epsilon, **sigma, **cut, **zeta, **mu, **beta;

  class AtomVecEllipsoid *avec;

  void allocate();
  double ylz_analytic(const int i, const int j, double a1[3][3], double a2[3][3], double *r12,
                      const double rsq, double *fforce, double *ttor, double *rtor);
};
}    // namespace LAMMPS_NS
#endif
#endif

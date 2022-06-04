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
PairStyle(brownian,PairBrownian);
// clang-format on
#else

#ifndef LMP_PAIR_BROWNIAN_H
#define LMP_PAIR_BROWNIAN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBrownian : public Pair {
 public:
  PairBrownian(class LAMMPS *);
  ~PairBrownian() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

 protected:
  double cut_inner_global, cut_global;
  double t_target, mu;
  int flaglog, flagfld;
  int flagHI, flagVF;
  int flagdeform, flagwall;
  double vol_P;
  double rad;
  class FixWall *wallfix;

  int seed;
  double **cut_inner, **cut;
  double R0, RT0;

  class RanMars *random;

  void set_3_orthogonal_vectors(double *, double *, double *);
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

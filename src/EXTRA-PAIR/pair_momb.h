/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
  -------------------------------------------------------------------------
  Contributed by Kristen Fichthorn, Tonnam Balankura, Ya Zhou @ Penn State University
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(momb,PairMomb);
// clang-format on
#else

#ifndef LMP_PAIR_MOMB_H
#define LMP_PAIR_MOMB_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMomb : public Pair {
 public:
  PairMomb(class LAMMPS *);
  ~PairMomb() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  double cut_global;
  double **cut;
  double sscale, dscale;
  double **d0, **alpha, **r0, **c, **rr;
  double **morse1;
  double **offset;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

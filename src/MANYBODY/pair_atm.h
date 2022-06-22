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
PairStyle(atm,PairATM);
// clang-format on
#else

#ifndef LMP_PAIR_ATM_H
#define LMP_PAIR_ATM_H

#include "pair.h"

namespace LAMMPS_NS {

class PairATM : public Pair {
 public:
  PairATM(class LAMMPS *);
  ~PairATM() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

 protected:
  double cut_global, cut_triple;
  double ***nu;

  void allocate();
  void interaction_ddd(double, double, double, double, double, double *, double *, double *,
                       double *, double *, int, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif

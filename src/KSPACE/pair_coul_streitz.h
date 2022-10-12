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
PairStyle(coul/streitz,PairCoulStreitz);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_Streitz_H
#define LMP_PAIR_COUL_Streitz_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulStreitz : public Pair {
 public:
  PairCoulStreitz(class LAMMPS *);
  ~PairCoulStreitz() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;
  void *extract(const char *, int &) override;

  static constexpr int NPARAMS_PER_LINE = 6;

 protected:
  struct Param {
    double chi, eta, gamma, zeta, zcore;
    int ielement;
  };

  int nmax;
  double precision;
  Param *params;    // parameter sets for elements

  // Kspace parameters
  int kspacetype;
  double cut_coul, cut_coulsq;
  double **scale;

  // Wolf
  double g_wolf, woself, dwoself;

  // Ewald
  double g_ewald;

  // QEq
  double *qeq_x, *qeq_j, *qeq_g, *qeq_z, *qeq_c;

  void allocate();
  virtual void read_file(char *);
  void setup_params();
  double self(Param *, double);
  void coulomb_integral_wolf(double, double, double, double &, double &, double &, double &);
  void wolf_sum(double, double, double, double, double, double, double, double, double &, double &);
  void coulomb_integral_ewald(double, double, double, double &, double &, double &, double &);
  void ewald_sum(double, double, double, double, double, double, double, double, double &, double &,
                 double);
};

}    // namespace LAMMPS_NS

#endif
#endif

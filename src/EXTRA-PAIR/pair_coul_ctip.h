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
PairStyle(coul/ctip,PairCoulCTIP);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_CTIP_H
#define LMP_PAIR_COUL_CTIP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulCTIP : public Pair {
 public:
  PairCoulCTIP(class LAMMPS *);
  ~PairCoulCTIP() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void *extract(const char *, int &) override;

  static constexpr int NPARAMS_PER_LINE = 9;

 protected:
  double cut_coul, cut_coulsq;
  double alpha;

  struct Param {
    double chi, eta, gamma, zeta, zcore, qmin, qmax, omega;
    int ielement;
  };

  Param *params;

  double **shield, **shieldcu, **reffc, **reffcsq, **reffc4, **reffc7;
  double **s2d_shift, **f_shift, **e_shift, **self_factor;

  double *qeq_x, *qeq_j, *qeq_g, *qeq_z, *qeq_c, *qeq_q1, *qeq_q2, *qeq_w;

  void allocate();

  void read_file(char *);
  void setup_params();
  double self(Param *, double);
};

}    // namespace LAMMPS_NS

#endif
#endif

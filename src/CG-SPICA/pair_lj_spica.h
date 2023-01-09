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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/spica,PairLJSPICA);
PairStyle(lj/sdk,PairLJSPICA);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_SPICA_H
#define LMP_PAIR_LJ_SPICA_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJSPICA : public Pair {
 public:
  PairLJSPICA(LAMMPS *);
  ~PairLJSPICA() override;
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
  void *extract(const char *, int &) override;
  double memory_usage() override;

 protected:
  int **lj_type;    // type of lennard jones potential

  double **cut;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset;

  // cutoff and offset for minimum of LJ potential
  // to be used in SPICA angle potential, which
  // uses only the repulsive part of the potential

  double **rminsq, **emin;

  double cut_global;

  virtual void allocate();

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();
};

}    // namespace LAMMPS_NS

#endif
#endif

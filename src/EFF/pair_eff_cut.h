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
PairStyle(eff/cut,PairEffCut);
// clang-format on
#else

#ifndef LMP_PAIR_EFF_CUT_H
#define LMP_PAIR_EFF_CUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairEffCut : public Pair {
 public:
  PairEffCut(class LAMMPS *);
  ~PairEffCut() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void min_pointers(double **, double **);
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  void min_xf_pointers(int, double **, double **) override;
  void min_xf_get(int) override;
  void min_x_set(int) override;
  double memory_usage() override;

 private:
  int limit_eradius_flag, pressure_with_evirials_flag;
  int ecp_found;
  double cut_global;
  double **cut;
  int ecp_type[100];
  double PAULI_CORE_A[100], PAULI_CORE_B[100], PAULI_CORE_C[100];
  double PAULI_CORE_D[100], PAULI_CORE_E[100];
  double hhmss2e, h2e;

  int nmax;
  double *min_eradius, *min_erforce;

  void allocate();
  void virial_eff_compute();
  void ev_tally_eff(int, int, int, int, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif

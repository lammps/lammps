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
PairStyle(lj/cut/coul/long/gauss,PairLJCutCoulLongGauss);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_LONG_GAUSS_H
#define LMP_PAIR_LJ_CUT_COUL_LONG_GAUSS_H

#include "electrode_pair.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutCoulLongGauss : public Pair, public ElectrodePair {

 public:
  PairLJCutCoulLongGauss(class LAMMPS *);
  ~PairLJCutCoulLongGauss() override;
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
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

  // ElectrodePair methods
  void compute_vector(double *, int, int, bool) override;
  void compute_vector_self(double *, int, int, bool) override;
  void compute_matrix(bigint *, double **, int) override;
  void compute_matrix_self(bigint *, double **, int) override;

 protected:
  double cut_lj_global;
  double **cut_lj, **cut_ljsq;
  double cut_coul, cut_coulsq;
  double **epsilon, **sigma, **eta;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  int *ispoint;
  double g_ewald;

  virtual void allocate();
 private:
  bool already_warned;
  void point_in_sensor_warning(int);
};

}    // namespace LAMMPS_NS

#endif
#endif

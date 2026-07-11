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
PairStyle(zbl,PairZBL);
// clang-format on
#else

#ifndef LMP_PAIR_ZBL_H
#define LMP_PAIR_ZBL_H

#include "pair.h"

namespace LAMMPS_NS {

class PairZBL : public Pair {
 public:
  PairZBL(class LAMMPS *);
  ~PairZBL() override;
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

 protected:
  double cut_global, cut_inner;
  double cut_globalsq, cut_innersq;
  double *z;
  double **d1a, **d2a, **d3a, **d4a, **zze;
  double **sw1, **sw2, **sw3, **sw4, **sw5;

  virtual void allocate();
  double e_zbl(double, int, int);
  double dzbldr(double, int, int);
  double d2zbldr2(double, int, int);
  void set_coeff(int, int, double, double);
};
}    // namespace LAMMPS_NS

#endif
#endif

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
PairStyle(table,PairTable);
// clang-format on
#else

#ifndef LMP_PAIR_TABLE_H
#define LMP_PAIR_TABLE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairTable : public Pair {
 public:
  PairTable(class LAMMPS *);
  ~PairTable() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

  enum { LOOKUP, LINEAR, SPLINE, BITMAP };

 protected:
  int tabstyle, tablength;
  struct Table {
    int ninput, rflag, fpflag, match, ntablebits;
    int nshiftbits, nmask;
    double rlo, rhi, fplo, fphi, cut;
    double *rfile, *efile, *ffile;
    double *e2file, *f2file;
    double innersq, delta, invdelta, deltasq6;
    double *rsq, *drsq, *e, *de, *f, *df, *e2, *f2;
  };
  int ntables;
  Table *tables;

  int **tabindex;

  virtual void allocate();
  void read_table(Table *, char *, char *);
  void param_extract(Table *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  virtual void compute_table(Table *);
  void null_table(Table *);
  void free_table(Table *);
  static void spline(double *, double *, int, double, double, double *);
  static double splint(double *, double *, double *, int, double);
};

}    // namespace LAMMPS_NS

#endif
#endif

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
PairStyle(multi/lucy,PairMultiLucy);
// clang-format on
#else

#ifndef LMP_PAIR_MULTI_LUCY_H
#define LMP_PAIR_MULTI_LUCY_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMultiLucy : public Pair {
 public:
  PairMultiLucy(class LAMMPS *);
  ~PairMultiLucy() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  void computeLocalDensity();
  double rho_0;

 protected:
  enum { LOOKUP, LINEAR };

  int nmax;

  int tabstyle, tablength;
  struct Table {
    int ninput, rflag, fpflag, match;
    double rlo, rhi, fplo, fphi, cut;
    double *rfile, *efile, *ffile;
    double *e2file, *f2file;
    double innersq, delta, invdelta, deltasq6;
    double *rsq, *drsq, *e, *de, *f, *df, *e2, *f2;
  };
  int ntables;
  Table *tables;

  int **tabindex;

  void allocate();
  void read_table(Table *, char *, char *);
  void param_extract(Table *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  void compute_table(Table *);
  void null_table(Table *);
  void free_table(Table *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);
};

}    // namespace LAMMPS_NS

#endif
#endif

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

#ifdef BOND_CLASS
// clang-format off
BondStyle(table,BondTable);
// clang-format on
#else

#ifndef LMP_BOND_TABLE_H
#define LMP_BOND_TABLE_H

#include "bond.h"

namespace LAMMPS_NS {

class BondTable : public Bond {
 public:
  BondTable(class LAMMPS *);
  ~BondTable() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  int tabstyle, tablength;
  double *r0;

  struct Table {
    int ninput, fpflag;
    double fplo, fphi, r0;
    double lo, hi;
    double *rfile, *efile, *ffile;
    double *e2file, *f2file;
    double delta, invdelta, deltasq6;
    double *r, *e, *de, *f, *df, *e2, *f2;
  };

  int ntables;
  Table *tables;
  int *tabindex;

  void allocate();
  void null_table(Table *);
  void free_table(Table *);
  void read_table(Table *, char *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  void compute_table(Table *);

  void param_extract(Table *, char *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);

  void uf_lookup(int, double, double &, double &);
  void u_lookup(int, double, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif

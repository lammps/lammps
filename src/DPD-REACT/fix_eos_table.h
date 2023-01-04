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

#ifdef FIX_CLASS
// clang-format off
FixStyle(eos/table,FixEOStable);
// clang-format on
#else

#ifndef LMP_FIX_EOS_TABLE_H
#define LMP_FIX_EOS_TABLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEOStable : public Fix {
 public:
  FixEOStable(class LAMMPS *, int, char **);
  ~FixEOStable() override;
  int setmask() override;
  void init() override;
  void post_integrate() override;
  void end_of_step() override;
  void energy_lookup(double, double &);
  void temperature_lookup(double, double &);

 protected:
  enum { LINEAR };

  int tabstyle, tablength;
  struct Table {
    int ninput;
    double lo, hi;
    double *rfile, *efile;
    double *e2file;
    double delta, invdelta, deltasq6;
    double *r, *e, *de, *e2;
  };
  int ntables;
  Table *tables;

  void allocate();
  void null_table(Table *);
  void free_table(Table *);
  void read_table(Table *, Table *, char *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  void compute_table(Table *);

  void param_extract(Table *, Table *, char *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);
};
}    // namespace LAMMPS_NS

#endif
#endif

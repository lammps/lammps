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
FixStyle(wall/table,FixWallTable);
// clang-format on
#else

#ifndef LMP_FIX_WALL_TABLE_H
#define LMP_FIX_WALL_TABLE_H

#include "fix_wall.h"

namespace LAMMPS_NS {

class FixWallTable : public FixWall {
 public:
  FixWallTable(class LAMMPS *, int, char **);
  ~FixWallTable() override;

  void post_constructor() override;
  void precompute(int) override;
  void wall_particle(int, int, double) override;

 protected:
  double offset[6];

  int tabstyle, tablength;
  struct Table {
    int ninput, fpflag;
    double fplo, fphi;
    double lo, hi;
    double *rfile, *efile, *ffile;
    double *e2file, *f2file;
    double delta, invdelta, deltasq6;
    double *r, *e, *de, *f, *df, *e2, *f2;
  };
  Table *tables;

  void null_table(Table &);
  void free_table(Table &);
  void read_table(Table &, const std::string &, const std::string &);
  void bcast_table(Table &);
  void spline_table(Table &);
  void compute_table(Table &);

  void param_extract(Table &, char *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);
  void uf_lookup(int, double, double &, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif

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
PairStyle(sw/angle/table,PairSWAngleTable);
// clang-format on
#else

#ifndef LMP_PAIR_SW_ANGLE_TABLE_H
#define LMP_PAIR_SW_ANGLE_TABLE_H

#include "pair_sw.h"

namespace LAMMPS_NS {

class PairSWAngleTable : public PairSW {
 public:
  PairSWAngleTable(class LAMMPS *);
  ~PairSWAngleTable() override;
  void compute(int, int) override;

  static constexpr int NPARAMS_PER_LINE = 18;

  // use struct Table as in class AngleTable
  struct Table {
    int ninput, fpflag;
    double fplo, fphi, theta0;
    double *afile, *efile, *ffile;
    double *e2file, *f2file;
    double delta, invdelta, deltasq6;
    double *ang, *e, *de, *f, *df, *e2, *f2;
  };

  struct ParamTable {
    int tablenamelength;        // length of table name
    char *tablename;            // name of associated angular table
    int keywordlength;          // length of key in table
    char *keyword;              // key in table
    int tabstyle, tablength;    // length of interpolation table (not ninput) and style
    Table *angtable;            // angle table
  };

 protected:
  ParamTable *table_params;    // tabulated parameter set for an I-J-K interaction

  void read_file(char *) override;
  void threebody_table(Param *, Param *, ParamTable *, double, double, double *, double *, double *,
                       double *, int, double &);

  void read_table(Table *, char *, char *);
  void spline_table(Table *);
  void compute_table(Table *, int length);
  void bcast_table(Table *);
  void null_table(Table *);
  void free_table(Table *);
  void free_param(ParamTable *);
  void param_extract(Table *, char *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);
  void uf_lookup(ParamTable *, double, double &, double &);
  void settings(int, char **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif

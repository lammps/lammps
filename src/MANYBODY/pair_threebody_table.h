/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(threebody/table,PairThreebodyTable);
// clang-format on
#else

#ifndef LMP_PAIR_THREEBODY_TABLE_H
#define LMP_PAIR_THREEBODY_TABLE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairThreebodyTable : public Pair {
 public:
  PairThreebodyTable(class LAMMPS *);
  ~PairThreebodyTable() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;

  static constexpr int NPARAMS_PER_LINE = 8;

  // no write or read from binary restart file

  // struct for threebody/table
  struct Table {
    int ninput;
    double rmin, rmax;
    double *r12file, *r13file, *thetafile, *f11file, *f12file, *f21file, *f22file, *f31file,
        *f32file, *efile;
  };

  struct Param {
    double cut, cutsq;
    int ielement, jelement, kelement;
    bool symmetric;             // whether it is a symmetric table or not
    int tablenamelength;        // length of table name
    char *tablename;            // name of associated angular table
    int keywordlength;          // length of key in table
    char *keyword;              // key in table
    int tabstyle, tablength;    // length of interpolation table (not ninput) and style
    Table *mltable;             // threebody table
  };

 protected:
  double cutmax;      // max cutoff for all elements
  Param *params;      // parameter set for an I-J-K interaction
  int maxshort;       // size of short neighbor list array
  int *neighshort;    // short neighbor list array

  void settings(int, char **) override;
  virtual void allocate();
  void read_file(char *);
  virtual void setup_params();
  void threebody(Param *, double, double, double *, double *, double *, double *, double *, int,
                 double &);

  void read_table(Table *, char *, char *, bool);
  void bcast_table(Table *, bool);
  void null_table(Table *);

  void free_table(Table *);
  void free_param(Param *);

  void param_extract(Table *, char *);

  void uf_lookup(Param *, double, double, double, double &, double &, double &, double &, double &,
                 double &, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif

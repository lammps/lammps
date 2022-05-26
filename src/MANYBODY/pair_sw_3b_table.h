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
PairStyle(sw/3b/table,PairSW3BTable);
// clang-format on
#else

#ifndef LMP_PAIR_SW_3B_TABLE_H
#define LMP_PAIR_SW_3B_TABLE_H

#include "pair_sw.h"

namespace LAMMPS_NS {

class PairSW3BTable : public PairSW {
 public:
  PairSW3BTable(class LAMMPS *);
  ~PairSW3BTable() override;
  void compute(int, int) override;

  static constexpr int NPARAMS_PER_LINE = 18;

  // use struct Table from class AngleTable
  struct Table {
    int ninput, fpflag;
    double fplo, fphi, theta0;
    double *afile, *efile, *ffile;
    double *e2file, *f2file;
    double delta, invdelta, deltasq6;
    double *ang, *e, *de, *f, *df, *e2, *f2;
  };
  
  struct Param_Extension {
    int tablenamelength;    // length of table name
    char *tablename;    // name of associated angular table
    int keywordlength;     // length of key in table
    char *keyword;    // key in table
    int tabstyle,tablength;    // length of interpolation table (not ninput) and style
    Table *angtable;    // angle table
  };
  
 protected:
  Param_Extension *extended_params;      // parameter set for an I-J-K interaction
   
  void read_file(char *) override;
  void threebody(Param *, Param *, Param *, Param_Extension*, double, double, double *, double *, double *,
                         double *, int, double &);

  void read_table(Table *, char *, char *);
  void spline_table(Table *);
  void compute_table(Table *,int length);
  void bcast_table(Table *);
  void null_table(Table *);
  void free_table(Table *);
  void free_param(Param_Extension *);
  void param_extract(Table *, char *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);
  void uf_lookup(Param_Extension *, double, double &, double &);                 
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style Stillinger-Weber requires atom IDs

This is a requirement to use the SW potential.

E: Pair style Stillinger-Weber requires newton pair on

See the newton command.  This is a restriction to use the SW
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open Stillinger-Weber potential file %s

The specified SW potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Stillinger-Weber potential file

Incorrect number of words per line in the potential file.

E: Illegal Stillinger-Weber parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

*/

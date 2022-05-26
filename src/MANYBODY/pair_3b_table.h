/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(3b/table,Pair3BTable)

#else

#ifndef LMP_PAIR_3B_TABLE_H
#define LMP_PAIR_3B_TABLE_H

#include "pair.h"

namespace LAMMPS_NS {

class Pair3BTable : public Pair {
 public:
  Pair3BTable(class LAMMPS *);
  ~Pair3BTable() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;

  static constexpr int NPARAMS_PER_LINE = 8;

  // no write or read from binary restart file

  // struct for 3b/table
  struct Table {
    int ninput;
    double rmin,rmax;
    double *r12file,*r13file,*thetafile,*f11file,*f12file,*f21file,*f22file,*f31file,*f32file,*efile;
  };

  struct Param {
    double cut,cutsq;
    int ielement,jelement,kelement;
    bool symmetric; // whether it is a symmetric table or not
    int tablenamelength; // length of table name
    char *tablename; // name of associated angular table
    int keywordlength; // length of key in table
    char *keyword; // key in table
    int tabstyle,tablength; // length of interpolation table (not ninput) and style
    Table *mltable; // 3b Table
  };

 protected:
  double cutmax;                // max cutoff for all elements
  Param *params;                // parameter set for an I-J-K interaction
  int maxshort;                 // size of short neighbor list array
  int *neighshort;              // short neighbor list array

  void settings(int, char **) override;
  virtual void allocate();
  void read_file(char *);
  virtual void setup_params();
  void threebody(Param *, Param *, Param *, double, double, double *, double *,
                 double *, double *, double *, int, double &);

  void read_table(Table *, char *, char *, bool);
  void bcast_table(Table *, bool);
  void null_table(Table *);

  void free_table(Table *);
  void free_param(Param *);

  void param_extract(Table *, char *);

  void uf_lookup(Param *, double, double, double, double &, double &, double &, double &, double &, double &, double &);
};

}


#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style 3b/table requires atom IDs

This is a requirement to use the 3b/table potential.

E: Pair style 3b/table requires newton pair on

See the newton command.  This is a restriction to use the 3b/table
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open 3b/table potential file %s

The specified 3b/table potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in 3b/table potential file

Incorrect number of words per line in the potential file.

E: Illegal 3b/table parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

*/

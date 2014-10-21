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

PairStyle(zbl,PairZBL)

#else

#ifndef LMP_PAIR_ZBL_H
#define LMP_PAIR_ZBL_H

#include "pair.h"

namespace LAMMPS_NS {

class PairZBL : public Pair {
 public:
  PairZBL(class LAMMPS *);
  virtual ~PairZBL();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_global, cut_inner;
  double cut_globalsq, cut_innersq;
  double *z;
  double **d1a,**d2a,**d3a,**d4a,**zze;
  double **sw1,**sw2,**sw3,**sw4,**sw5;

  void allocate();
  double e_zbl(double, int, int);
  double dzbldr(double, int, int);
  double d2zbldr2(double, int, int);
};

namespace PairZBLConstants {

  // ZBL constants

  static const double pzbl = 0.23;
  static const double a0 = 0.46850;
  static const double c1 = 0.02817;
  static const double c2 = 0.28022;
  static const double c3 = 0.50986;
  static const double c4 = 0.18175;
  static const double d1 = 0.20162;
  static const double d2 = 0.40290;
  static const double d3 = 0.94229;
  static const double d4 = 3.19980;
}

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

*/

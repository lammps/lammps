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
PairStyle(zbl,PairZBL);
// clang-format on
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
  virtual void init_style();
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_global, cut_inner;
  double cut_globalsq, cut_innersq;
  double *z;
  double **d1a, **d2a, **d3a, **d4a, **zze;
  double **sw1, **sw2, **sw3, **sw4, **sw5;

  virtual void allocate();
  double e_zbl(double, int, int);
  double dzbldr(double, int, int);
  double d2zbldr2(double, int, int);
  void set_coeff(int, int, double, double);
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

*/

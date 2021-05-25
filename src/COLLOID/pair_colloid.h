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
PairStyle(colloid,PairColloid);
// clang-format on
#else

#ifndef LMP_PAIR_COLLOID_H
#define LMP_PAIR_COLLOID_H

#include "pair.h"

namespace LAMMPS_NS {

class PairColloid : public Pair {
 public:
  PairColloid(class LAMMPS *);
  virtual ~PairColloid();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 protected:
  enum { SMALL_SMALL, SMALL_LARGE, LARGE_LARGE };

  double cut_global;
  double **cut;
  double **a12, **d1, **d2, **diameter, **a1, **a2, **offset;
  double **sigma, **sigma3, **sigma6;
  double **lj1, **lj2, **lj3, **lj4;
  int **form;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Overlapping small/large in pair colloid

This potential is infinite when there is an overlap.

E: Overlapping large/large in pair colloid

This potential is infinite when there is an overlap.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Invalid d1 or d2 value for pair colloid coeff

Neither d1 or d2 can be < 0.

*/

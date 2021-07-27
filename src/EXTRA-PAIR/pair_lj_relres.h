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
PairStyle(lj/relres,PairLJRelRes);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_RELRES_H
#define LMP_PAIR_LJ_RELRES_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJRelRes : public Pair {
 public:
  PairLJRelRes(class LAMMPS *);
  virtual ~PairLJRelRes();
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
  double cut_inner_global, cut_global, cutf_inner_global, cutf_global;
  double **cut, **cut_inner, **cut_inner_sq, **cutf, **cutfsq, **cutf_inner, **cutf_inner_sq;
  double **epsilon, **sigma;
  double **epsilonf, **sigmaf;
  double **lj1, **lj2, **lj3, **lj4;
  double **ljf1, **ljf2, **ljf3, **ljf4;
  double **ljsw0, **ljsw1, **ljsw2, **ljsw3, **ljsw4;
  double **ljswf0, **ljswf1, **ljswf2, **ljswf3, **ljswf4;
  double **offset, **offsetsp, **offsetsm;

  void allocate();
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

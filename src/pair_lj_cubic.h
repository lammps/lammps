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

PairStyle(lj/cubic,PairLJCubic)

#else

#ifndef LMP_PAIR_LJ_CUBIC_H
#define LMP_PAIR_LJ_CUBIC_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCubic : public Pair {
 public:
  PairLJCubic(class LAMMPS *);
  virtual ~PairLJCubic();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 protected:
  double **cut,**cut_inner,**cut_inner_sq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4;

  void allocate();
};

namespace PairLJCubicConstants {

  // LJ quantities scaled by epsilon and rmin = sigma*2^1/6

  static const double RT6TWO = 1.1224621;  // 2^1/6
  static const double SS = 1.1086834;      // inflection point (13/7)^1/6
  static const double PHIS = -0.7869823;   // energy at s
  static const double DPHIDS = 2.6899009;  // gradient at s
  static const double A3 = 27.93357;       // cubic coefficient
  static const double SM = 1.5475375;      // cubic cutoff = s*67/48

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

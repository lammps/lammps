/* ----------------------------------------------------------------------
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
  virtual void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double **cut,**cut_inner,**cut_inner_sq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4;

  // LJ quantities scaled by epsilon and rmin = sigma*2^1/6

  static const double rt6two = 1.1224621;  // 2^1/6
  static const double s = 1.1086834;       // inflection point = (13/7)^1/6
  static const double phis = -0.7869823;   // energy at s
  static const double dphids = 2.6899009;  // gradient at s
  static const double a3 = 27.93357;       // cubic coefficient
  static const double sm = 1.5475375;      // cubic cutoff = s*67/48

  void allocate();
};

}

#endif
#endif

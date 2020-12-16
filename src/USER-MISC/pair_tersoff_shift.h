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

PairStyle(tersoff/shift,PairTersoffShift)

#else

#ifndef LMP_PAIR_TERSOFF_SHIFT_H
#define LMP_PAIR_TERSOFF_SHFIT_H

#include "pair_tersoff.h"

namespace LAMMPS_NS {

class PairTersoffShift : public PairTersoff {
 public:
  PairTersoffShift(class LAMMPS *);
  ~PairTersoffShift() {}
//  void settings(int, char **);

// protected:
//  virtual void repulsive(Param *, double, double &, int, double &);
//  virtual double zeta(Param *, double, double, double *, double *);
//  virtual void force_zeta(Param *, double, double, double &,
//                          double &, int, double &);
//  void attractive(Param *, double, double, double, double *, double *,
//                  double *, double *, double *);
//  void costheta_d(double *, double, double *, double,
//                  double *, double *, double *);
//
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair tersoff/zbl requires metal or real units

This is a current restriction of this pair potential.

E: Cannot open Tersoff potential file %s

The specified potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Tersoff potential file

Incorrect number of words per line in the potential file.

E: Illegal Tersoff parameter

One or more of the coefficients defined in the potential file is
invalid.

*/

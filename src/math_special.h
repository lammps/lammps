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

#ifndef LMP_MATH_SPECIAL_H
#define LMP_MATH_SPECIAL_H

#include <math.h>

namespace LAMMPS_NS {

namespace MathSpecial {
  // x**2;
  static inline double square(const double &x) { return x*x; }

  // optimized version of (sin(x)/x)**n with n being a positive integer
  static inline double powsinxx(const double x, int n) {
    double xx,yy,ww;

    if ((x == 0.0) || (n == 0)) return 1.0;

    xx = sin(x)/x;
    yy = (n & 1) ? xx : 1.0;
    ww = xx;
    n >>= 1;

    while (n) {
      ww *= ww;
      if (n & 1) yy *= ww;
      n >>=1;
    }
    return yy;
  }
}
}

#endif

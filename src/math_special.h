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

  // x**2, use instead of pow(x,2.0)

  static inline double square(const double &x) { return x*x; }

  // x**3, use instead of pow(x,3.0)
  static inline double cube(const double &x) { return x*x*x; }

  // return -1.0 for odd n, 1.0 for even n, like pow(-1.0,n)
  static inline double powsign(const int n) { return (n & 1) ? -1.0 : 1.0; }

  // optimized version of pow(x,n) with n being integer
  // up to 10x faster than pow(x,y)

  static inline double powint(const double &x, const int n) {
    double yy,ww;

    if (x == 0.0) return 0.0;
    int nn = (n > 0) ? n : -n;
    ww = x;

    for (yy = 1.0; nn != 0; nn >>= 1, ww *=ww)
      if (nn & 1) yy *= ww;

    return (n > 0) ? yy : 1.0/yy;
  }

  // optimized version of (sin(x)/x)**n with n being a _positive_ integer

  static inline double powsinxx(const double &x, int n) {
    double yy,ww;

    if (x == 0.0) return 1.0;

    ww = sin(x)/x;

    for (yy = 1.0; n != 0; n >>= 1, ww *=ww)
      if (n & 1) yy *= ww;

    return yy;
  }
}
}

#endif

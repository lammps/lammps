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

/* ----------------------------------------------------------------------
   Contributing author: Tim Bernhard (ETHZ)
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(fourier/simple/approx,AngleFourierSimpleApprox);
// clang-format on
#else

#ifndef LMP_ANGLE_FOURIER_SIMPLE_APPROX_H
#define LMP_ANGLE_FOURIER_SIMPLE_APPROX_H

#include "angle_fourier_simple.h"
#include "math_const.h"

namespace LAMMPS_NS {

   /* ---------------------------------------------------------------------- */

/**
 * Approximation of cos
 * maximum error is about 0.00109 for the range -pi to pi
 * Source: https://stackoverflow.com/a/28050328/3909202
 */
static double fastCos(double x) {
  constexpr double inv2pi = 1./(MathConst::MY_2PI);
  // TODO: check if range map is necessary
  double x_wrapped = x - MathConst::MY_2PI * floor(x * inv2pi);
    x_wrapped *= inv2pi;
    x_wrapped -= 0.25 + floor(x_wrapped + 0.25);
    x_wrapped *= 16.0 * (abs(x_wrapped) - 0.5);
    x_wrapped += 0.225 * x_wrapped * (abs(x_wrapped) - 1.0);
    return x_wrapped;
}

/* ---------------------------------------------------------------------- */

/**
 * Approximation of acos 
 * Returns the arccosine of x in the range [0,pi], expecting x to be in the range [-1,+1]. 
 * Absolute error <= 6.7e-5
 * Range mapping not necessary as C++ standard library would throw domain error
 * Source: https://developer.download.nvidia.com/cg/acos.html
 */
static double fastAcos(double x) {
  double negate = double(x < 0);
  x = fabs(x);
  double ret = -0.0187293;
  ret = ret * x;
  // TODO: verify that fpas are used here
  ret = ret + 0.0742610;
  ret = ret * x;
  ret = ret - 0.2121144;
  ret = ret * x;
  ret = ret + 1.5707288;
  ret = ret * sqrt(1.0-x); 
  ret = ret - 2 * negate * ret;
  return negate * 3.14159265358979 + ret;
}

// static double fastCos(double);
// static double fastAcos(double);

class AngleFourierSimpleApprox : public AngleFourierSimple {

 public:
  AngleFourierSimpleApprox(class LAMMPS *lmp);
  virtual void compute(int, int);

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND>
  void eval();
};

}    // namespace LAMMPS_NS

#endif
#endif

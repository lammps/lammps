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
   [ based on angle_fourier_simple_omp.cpp Axel Kohlmeyer (Temple U)]
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
#include <cassert>

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */

#define ONE_OVER_PI 0.318309886183790671538035357
#define ONE_OVER_PI_SQUARE 0.1013211836423377714440499
#define ONE_OVER_PI_CUBED 0.03225153443319948918450346
#define TWO_PI_INVERSE 0.1591549430918953357690176

/**
 * Approximation of sine. |x*e(x)|≤2e-9
 * @source M. Abramowitz and I. A. Stegun, Eds., Handbook of mathematical functions: with formulas, graphs, and mathematical tables, p. 43
 */
static float sinePolynomial(float x)
{
  const float x2 = x * x;

  const float term1 = x2 * -0.0000000239 + 0.0000027526;
  const float term2 = x2 * term1 - 0.0001984090;
  const float term3 = x2 * term2 + 0.00833333315;
  const float term4 = x2 * term3 - 0.16666666664;
  return x * (1 + x2 * term4);
}

/**
 * Approximation of sine
 * @source M. Abramowitz and I. A. Stegun, Eds., Handbook of mathematical functions: with formulas, graphs, and mathematical tables, p. 43
 */
static float fastSin(float x)
{
  // wrap x within [0, TWO_PI)
  const float a = x * TWO_PI_INVERSE;
  x -= static_cast<int>(a) * MathConst::MY_2PI;
  if (x < 0.0f) x += MathConst::MY_2PI;
  // assert(x >= 0);
  // assert(x < MathConst::MY_2PI);

  // 4 pieces of hills: wrap x within [0, pi/2]
  if (x < MathConst::MY_PI2)
    return sinePolynomial(x);
  else if (x < MathConst::MY_PI)
    return sinePolynomial(MathConst::MY_PI - x);
  else if (x < 3.0f * MathConst::MY_PI2)
    return -sinePolynomial(x - MathConst::MY_PI);
  else
    return -sinePolynomial(MathConst::MY_2PI - x);
}

/**
 * Approximation of cosine. |e(x)|≤2e-9
 * @source M. Abramowitz and I. A. Stegun, Eds., Handbook of mathematical functions: with formulas, graphs, and mathematical tables, p. 43
 */
static float cosinePolynomial(float x)
{
  const float x2 = x * x;

  const float term1 = x2 * -0.0000002605 + 0.0000247609;
  const float term2 = x2 * term1 - 0.0013888397;
  const float term3 = x2 * term2 + 0.0416666418;
  const float term4 = x2 * term3 - 0.4999999963;
  return 1.0 + x2 * term4;
}

/**
 * Approximation of cosine
 * @source M. Abramowitz and I. A. Stegun, Eds., Handbook of mathematical functions: with formulas, graphs, and mathematical tables, p. 43
 */
static float fastCos(float x)
{
  // wrap x within [0, TWO_PI)
  const float a = x * TWO_PI_INVERSE;
  x -= static_cast<int>(a) * MathConst::MY_2PI;
  if (x < 0.0f) x += MathConst::MY_2PI;
  // const float x = fmodf(MathConst::MY_2PI + fmodf(a, MathConst::MY_2PI), MathConst::MY_2PI);

  // 4 pieces of hills: wrap x within [0, pi/2]
  if (x < MathConst::MY_PI2)
    return cosinePolynomial(x);
  else if (x < MathConst::MY_PI)
    return -cosinePolynomial(MathConst::MY_PI - x);
  else if (x < 3.0f * MathConst::MY_PI2)
    return -cosinePolynomial(x - MathConst::MY_PI);
  else
    return cosinePolynomial(MathConst::MY_2PI - x);
}

/* ---------------------------------------------------------------------- */

/**
 * Combined call to sin & cos. Gets optimized by compiler where possibly to simultaneous call.
 * More stable/portable than <math.h>'s call, it seems.
 */
std::pair<double,double> sincos(double arg) { return { std::sin(arg), std::cos(arg) }; }

/* ---------------------------------------------------------------------- */

/**
 * Approximation of acos 
 * Returns the arccosine of x in the range [0,pi], expecting x to be in the range [-1,+1]. 
 * Absolute error <= 6.7e-5
 * Range mapping not necessary as C++ standard library would throw domain error
 * Source: https://developer.download.nvidia.com/cg/acos.html
 */
static double fastAcos(double x)
{
  double negate = double(x < 0);
  x = fabs(x);
  // fmas are used here only if I use -ffast-math
  // tested compiler: clang version 12.0.0, Target: x86_64-apple-darwin20.5.0
  double ret = -0.0187293 * x + 0.0742610;
  double ret2 = ret * x - 0.2121144;
  double ret3 = ret2 * x + 1.5707288;
  double ret4 = ret3 * sqrt(1.0 - x);
  double ret5 = ret4 * (1 - 2 * negate);
  return negate * 3.14159265358979 + ret5;
}

// static double fastCos(double);
// static double fastAcos(double);

class AngleFourierSimpleApprox : public AngleFourierSimple {

 public:
  AngleFourierSimpleApprox(class LAMMPS *lmp);
  virtual void compute(int, int);

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND> void eval();
};

}    // namespace LAMMPS_NS

#endif
#endif

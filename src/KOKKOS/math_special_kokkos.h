// clang-format off
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

#ifndef LMP_MATH_SPECIAL_KOKKOS_H
#define LMP_MATH_SPECIAL_KOKKOS_H

#include <cmath>
#include "kokkos_type.h"

namespace LAMMPS_NS {

namespace MathSpecialKokkos {

  /*! Fast tabulated factorial function
   *
   *  This function looks up pre-computed factorial values for arguments of n = 0
   *  to a maximum of 167, which is the maximal value representable by a double
   *  precision floating point number.  For other values of n a NaN value is returned.
   *
   *  \param   n  argument (valid: 0 <= n <= 167)
   *  \return  value of n! as double precision number or NaN */

  extern double factorial(const int n);

  /* optimizer friendly implementation of exp2(x).
   *
   * strategy:
   *
   * split argument into an integer part and a fraction:
   * ipart = floor(x+0.5);
   * fpart = x - ipart;
   *
   * compute exp2(ipart) from setting the ieee754 exponent
   * compute exp2(fpart) using a pade' approximation for x in [-0.5;0.5[
   *
   * the result becomes: exp2(x) = exp2(ipart) * exp2(fpart)
   */

  /* IEEE 754 double precision floating point data manipulation */
  typedef union
  {
    double   f;
    uint64_t u;
    struct {int32_t  i0,i1;} s;
  }  udi_t;

  /* double precision constants */
  #define FM_DOUBLE_LOG2OFE  1.4426950408889634074

  /*! Fast implementation of 2^x without argument checks for little endian CPUs
   *
   *  This function implements an optimized version of pow(2.0, x) that does not
   *  check for valid arguments and thus may only be used where arguments are well
   *  behaved.  The implementation makes assumptions about the layout of double
   *  precision floating point numbers in memory and thus will only work on little
   *  endian CPUs.  If little endian cannot be safely detected, the result of
   *  calling pow(2.0, x) will be returned.  This function also is the basis for
   *  the fast exponential fm_exp(x).
   *
   *  \param   x argument
   *  \return  value of 2^x as double precision number */

  KOKKOS_INLINE_FUNCTION
  static double exp2_x86(double x)
  {
  #if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
      double   ipart, fpart, px, qx;
      udi_t    epart;

  const double fm_exp2_q[2] = {
  /*  1.00000000000000000000e0, */
      2.33184211722314911771e2,
      4.36821166879210612817e3
  };
  const double fm_exp2_p[3] = {
      2.30933477057345225087e-2,
      2.02020656693165307700e1,
      1.51390680115615096133e3
  };

      ipart = floor(x+0.5);
      fpart = x - ipart;
      epart.s.i0 = 0;
      epart.s.i1 = (((int) ipart) + 1023) << 20;

      x = fpart*fpart;

      px =        fm_exp2_p[0];
      px = px*x + fm_exp2_p[1];
      qx =    x + fm_exp2_q[0];
      px = px*x + fm_exp2_p[2];
      qx = qx*x + fm_exp2_q[1];

      px = px * fpart;

      x = 1.0 + 2.0*(px/(qx-px));
      return epart.f*x;
  #else
      return pow(2.0, x);
  #endif
  }

  /*! Fast implementation of exp(x) for little endian CPUs
   *
   *  This function implements an optimized version of exp(x) for little endian CPUs.
   *  It calls the exp2_x86(x) function with a suitable prefactor to x to return exp(x).
   *  The implementation makes assumptions about the layout of double
   *  precision floating point numbers in memory and thus will only work on little
   *  endian CPUs.  If little endian cannot be safely detected, the result of
   *  calling the exp(x) implementation in the standard math library will be returned.
   *
   *  \param   x argument
   *  \return  value of e^x as double precision number */

  KOKKOS_INLINE_FUNCTION
  static double fm_exp(double x)
  {
  #if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
      if (x < -1022.0/FM_DOUBLE_LOG2OFE) return 0;
      if (x > 1023.0/FM_DOUBLE_LOG2OFE) return INFINITY;
      return exp2_x86(FM_DOUBLE_LOG2OFE * x);
  #else
      return ::exp(x);
  #endif
  }

  // support function for scaled error function complement

  extern double erfcx_y100(const double y100);

  /*! Fast scaled error function complement exp(x*x)*erfc(x) for coul/long styles
   *
   *  This is a portable fast implementation of exp(x*x)*erfc(x) that can be used
   *  in coul/long pair styles as a replacement for the polynomial expansion that
   *  is/was widely used.  Unlike the polynomial expansion, that is only accurate
   *  at the level of single precision floating point it provides full double precision
   *  accuracy, but at comparable speed (unlike the erfc() implementation shipped
   *  with GNU standard math library).
   *
   *  \param   x argument
   *  \return  value of e^(x*x)*erfc(x) */

  static inline double my_erfcx(const double x)
  {
    if (x >= 0.0)
      return erfcx_y100(400.0 / (4.0 + x));
    else
      return 2.0 * exp(x * x) - erfcx_y100(400.0 / (4.0 - x));
  }

  /*! Fast implementation of exp(-x*x) for little endian CPUs for coul/long styles
   *
   *  This function implements an optimized version of exp(-x*x) based on exp2_x86()
   *  for use with little endian CPUs. If little endian cannot be safely detected,
   *  the result of calling the exp(-x*x) implementation in the standard math
   *  library will be returned.
   *
   *  \param   x argument
   *  \return  value of e^(-x*x) as double precision number */

  static inline double expmsq(double x)
  {
    x *= x;
    x *= 1.4426950408889634074; // log_2(e)
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return (x < 1023.0) ? exp2_x86(-x) : 0.0;
#else
    return (x < 1023.0) ? exp2(-x) : 0.0;
#endif
  }

  /*! Fast inline version of pow(x, 2.0)
   *
   *  \param   x argument
   *  \return  x*x */

  KOKKOS_INLINE_FUNCTION
  static double square(const double &x) { return x * x; }

  /*! Fast inline version of pow(x, 3.0)
   *
   *  \param   x argument
   *  \return  x*x */

  KOKKOS_INLINE_FUNCTION
  static double cube(const double &x) { return x * x * x; }

  /* Fast inline version of pow(-1.0, n)
   *
   *  \param   n argument (integer)
   *  \return  -1 if n is odd, 1.0 if n is even */

  KOKKOS_INLINE_FUNCTION
  static double powsign(const int n) { return (n & 1) ? -1.0 : 1.0; }

  /* Fast inline version of pow(x,n) for integer n
   *
   * This is a version of pow(x,n) optimized for n being integer.
   * Speedups of up to 10x faster than pow(x,y) have been measured.
   *
   *  \param   n argument (integer)
   *  \return  value of x^n */

  KOKKOS_INLINE_FUNCTION
  static double powint(const double &x, const int n)
  {
    double yy, ww;

    if (x == 0.0) return 0.0;
    int nn = (n > 0) ? n : -n;
    ww = x;

    for (yy = 1.0; nn != 0; nn >>= 1, ww *= ww)
      if (nn & 1) yy *= ww;

    return (n > 0) ? yy : 1.0 / yy;
  }

  /* Fast inline version of (sin(x)/x)^n as used by PPPM kspace styles
   *
   * This is an optimized function to compute (sin(x)/x)^n as frequently used by PPPM.
   *
   *  \param   n argument (integer). Expected to be positive.
   *  \return  value of (sin(x)/x)^n */

  KOKKOS_INLINE_FUNCTION
  static double powsinxx(const double &x, int n)
  {
    double yy, ww;

    if (x == 0.0) return 1.0;

    ww = sin(x) / x;

    for (yy = 1.0; n != 0; n >>= 1, ww *= ww)
      if (n & 1) yy *= ww;

    return yy;
  }
}    // namespace MathSpecialKokkos
}    // namespace LAMMPS_NS

#endif

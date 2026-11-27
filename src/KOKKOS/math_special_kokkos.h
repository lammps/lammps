// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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



namespace LAMMPS_NS::MathSpecialKokkos {

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

  /*! Fast implementation of exp(-x*x) for coul/long styles
   *
   *  This function implements an optimized version of exp(-x*x)
   *
   *  \param   x argument
   *  \return  value of e^(-x*x) as double precision number */

  static inline double expmsq(double x)
  {
    x *= x;
    return (x < 1023.0) ? Kokkos::exp(-x) : 0.0;
  }

  /*! Fast inline version of pow(x, 2.0)
   *
   *  \param   x argument
   *  \return  x*x */
  template<typename T>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static T square(const T &x) { return x * x; }

  /*! Fast inline version of pow(x, 3.0)
   *
   *  \param   x argument
   *  \return  x*x */
  template<typename T>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static T cube(const T &x) { return x * x * x; }

  /* Fast inline version of pow(-1.0, n)
   *
   *  \param   n argument (integer)
   *  \return  -1 if n is odd, 1.0 if n is even */

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static KK_FLOAT powsign(const int n) { return (n & 1) ? static_cast<KK_FLOAT>(-1.0) : static_cast<KK_FLOAT>(1.0); }

  /* Fast inline version of pow(x,n) for integer n
   *
   * This is a version of pow(x,n) optimized for n being integer.
   * Speedups of up to 10x faster than pow(x,y) have been measured.
   *
   *  \param   n argument (integer)
   *  \return  value of x^n */

  template<typename T>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static T powint(const T &x, const int n)
  {
    T yy, ww;

    if (x == static_cast<T>(0)) return static_cast<T>(0);
    int nn = (n > 0) ? n : -n;
    ww = x;

    for (yy = static_cast<T>(1); nn != 0; nn >>= 1, ww *= ww)
      if (nn & 1) yy *= ww;

    return (n > 0) ? yy : static_cast<T>(1) / yy;
  }

  /* Version of pow(x,n) for an exponent that is only known at run time
   *
   * Device capable equivalent of MathSpecial::powauto().  Dispatches to
   * powint() when the exponent happens to have an integer value within the
   * range of an int, and to Kokkos::pow() otherwise.  Follows the powint()
   * convention of returning 0 for x == 0, which differs from pow() for a
   * negative exponent.
   *
   *  \param   x base
   *  \param   n exponent
   *  \return  value of x^n */

  template<typename T>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static T powauto(const T &x, const T &n)
  {
    if (n == static_cast<T>(0)) return static_cast<T>(1);
    if (x == static_cast<T>(0)) return static_cast<T>(0);
    const T rounded = Kokkos::round(n);
    // the conversion to int is only defined for |n| < 2^31, which excludes inf
    if ((n == rounded) && (Kokkos::fabs(rounded) < static_cast<T>(2147483648.0)))
      return powint(x, static_cast<int>(rounded));
    return Kokkos::pow(x, n);
  }

  /* Fast inline version of (sin(x)/x)^n as used by PPPM kspace styles
   *
   * This is an optimized function to compute (sin(x)/x)^n as frequently used by PPPM.
   *
   *  \param   n argument (integer). Expected to be positive.
   *  \return  value of (sin(x)/x)^n */

  template<typename T>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static T powsinxx(const T &x, int n)
  {
    T yy, ww;

    if (x == static_cast<T>(0)) return static_cast<T>(1);

    ww = sin(x) / x;

    for (yy = static_cast<T>(1); n != 0; n >>= 1, ww *= ww)
      if (n & 1) yy *= ww;

    return yy;
  }

  /* ----------------------------------------------------------------------
    ans = v1 - v2
  ------------------------------------------------------------------------- */

  template<typename T>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static void sub3(const T *v1, const T *v2, T *ans)
  {
    ans[0] = v1[0] - v2[0];
    ans[1] = v1[1] - v2[1];
    ans[2] = v1[2] - v2[2];
  }

  /* ----------------------------------------------------------------------
    dot product of 2 vectors
  ------------------------------------------------------------------------- */

  template<typename T>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static T dot3(const T *v1, const T *v2)
  {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  }

  /* ----------------------------------------------------------------------
    ans = v1 cross v2
  ------------------------------------------------------------------------- */

  template<typename T>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static void cross3(const T *v1, const T *v2, T *ans)
  {
    ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
    ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
    ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }

} // namespace LAMMPS_NS::MathSpecialKokkos


#endif

//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_MATHEMATICAL_SPECIAL_FUNCTIONS_HPP
#define KOKKOS_MATHEMATICAL_SPECIAL_FUNCTIONS_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_MATHSPECFUNCTIONS
#endif

#include <Kokkos_Macros.hpp>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <Kokkos_NumericTraits.hpp>
#include <Kokkos_Complex.hpp>

namespace Kokkos {
namespace Experimental {

//! Compute exponential integral E1(x) (x > 0).
template <class RealType>
KOKKOS_INLINE_FUNCTION RealType expint1(RealType x) {
  // This function is a conversion of the corresponding Fortran program in
  // S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
  using Kokkos::exp;
  using Kokkos::fabs;
  using Kokkos::log;
  using Kokkos::pow;
  using Kokkos::Experimental::epsilon;
  using Kokkos::Experimental::infinity;

  RealType e1;

  if (x < 0) {
    e1 = -infinity<RealType>::value;
  } else if (x == 0.0) {
    e1 = infinity<RealType>::value;
  } else if (x <= 1.0) {
    e1         = 1.0;
    RealType r = 1.0;
    for (int k = 1; k <= 25; k++) {
      RealType k_real = static_cast<RealType>(k);
      r               = -r * k_real * x / pow(k_real + 1.0, 2.0);
      e1              = e1 + r;
      if (fabs(r) <= fabs(e1) * epsilon<RealType>::value) break;
    }
    e1 = -0.5772156649015328 - log(x) + x * e1;
  } else {
    int m       = 20 + static_cast<int>(80.0 / x);
    RealType t0 = 0.0;
    for (int k = m; k >= 1; k--) {
      RealType k_real = static_cast<RealType>(k);
      t0              = k_real / (1.0 + k_real / (x + t0));
    }
    e1 = exp(-x) * (1.0 / (x + t0));
  }
  return e1;
}

//! Compute error function erf(z) for z=cmplx(x,y).
template <class RealType>
KOKKOS_INLINE_FUNCTION Kokkos::complex<RealType> erf(
    const Kokkos::complex<RealType>& z) {
  // This function is a conversion of the corresponding Fortran program written
  // by D.E. Amos, May,1974. D.E. Amos' revisions of Jan 86 incorporated by
  // Ken Damrau on 27-Jan-1986 14:37:13
  //
  // Reference: NBS HANDBOOK OF MATHEMATICAL FUNCTIONS, AMS 55, By
  //           M. ABRAMOWITZ AND I.A. STEGUN, December,1955.
  // Summary:
  //  If x < 0, z is replaced by -z and all computation is done in the right
  //  half lane, except for z inside the circle abs(z)<=2, since
  //  erf(-z)=-erf(z). The regions for computation are divided as follows
  //      (1)  abs(z)<=2 - Power series, NBS Handbook, p. 298
  //      (2)  abs(z)>2 and x>1 - continued fraction, NBS Handbook, p. 298
  //      (3)  abs(z)>2 and 0<=x<=1 and abs(y)<6 - series, NBS Handbook, p. 299
  //      (4)  abs(z)>2 and 0<=x<=1 and abs(y)>=6 - asymptotic expansion
  //  Error condition: abs(z^2) > 670 is a fatal overflow error
  using Kokkos::cos;
  using Kokkos::exp;
  using Kokkos::fabs;
  using Kokkos::sin;
  using Kokkos::Experimental::epsilon_v;
  using Kokkos::Experimental::infinity_v;
  using Kokkos::numbers::pi_v;

  using CmplxType = Kokkos::complex<RealType>;

  constexpr auto inf = infinity_v<RealType>;
  constexpr auto tol = epsilon_v<RealType>;

  const RealType fnorm = 1.12837916709551;
  const RealType gnorm = 0.564189583547756;
  const RealType eh    = 0.606530659712633;
  const RealType ef    = 0.778800783071405;
  // const RealType tol   = 1.0e-13;
  constexpr auto pi = pi_v<RealType>;

  CmplxType cans;

  RealType az = Kokkos::abs(z);
  if (az <= 2.0) {  // Series for abs(z)<=2.0
    CmplxType cz    = z * z;
    CmplxType accum = CmplxType(1.0, 0.0);
    CmplxType term  = accum;
    RealType ak     = 1.5;
    for (int i = 1; i <= 35; i++) {
      term  = term * cz / ak;
      accum = accum + term;
      if (Kokkos::abs(term) <= tol) break;
      ak = ak + 1.0;
    }
    cz          = -cz;
    RealType er = cz.real();
    RealType ei = cz.imag();
    accum       = accum * z * fnorm;
    cz          = exp(er) * CmplxType(cos(ei), sin(ei));
    cans        = accum * cz;
  }       // end (az <= 2.0)
  else {  //(az > 2.0)
    CmplxType zp = z;
    if (z.real() < 0.0) zp = -z;
    CmplxType cz = zp * zp;
    RealType xp  = zp.real();
    RealType yp  = zp.imag();
    if (xp > 1.0) {
      // continued fraction for erfc(z), abs(Z)>2
      int n          = static_cast<int>(100.0 / az + 5.0);
      int fn         = n;
      CmplxType term = cz;
      for (int i = 1; i <= n; i++) {
        RealType fnh = fn - 0.5;
        term         = cz + (fnh * term) / (fn + term);
        fn           = fn - 1;
      }
      if (Kokkos::abs(cz) > 670.0) return CmplxType(inf, inf);
      cz              = -cz;
      RealType er     = cz.real();
      RealType ei     = cz.imag();
      cz              = exp(er) * CmplxType(cos(ei), sin(ei));
      CmplxType accum = zp * gnorm * cz;
      cans            = 1.0 - accum / term;
      if (z.real() < 0.0) cans = -cans;
    }       // end (xp > 1.0)
    else {  //(xp <= 1.0)
      if (fabs(yp) <
          6.0) {  // Series (3) for abs(z)>2 and 0<=xp<=1 and abs(yp)<6
        RealType s1   = 0.0;
        RealType s2   = 0.0;
        RealType x2   = xp * xp;
        RealType fx2  = 4.0 * x2;
        RealType tx   = xp + xp;
        RealType xy   = xp * yp;
        RealType sxyh = sin(xy);
        RealType sxy  = sin(xy + xy);
        RealType cxy  = cos(xy + xy);
        RealType fn   = 1.0;
        RealType fnh  = 0.5;
        RealType ey   = exp(yp);
        RealType en   = ey;
        RealType ehn  = eh;
        RealType un   = ef;
        RealType vn   = 1.0;
        for (int i = 1; i <= 50; i++) {
          RealType ren = 1.0 / en;
          RealType csh = en + ren;
          RealType tm  = xp * csh;
          RealType ssh = en - ren;
          RealType tmp = fnh * ssh;
          RealType rn  = tx - tm * cxy + tmp * sxy;
          RealType ain = tm * sxy + tmp * cxy;
          RealType cf  = un / (vn + fx2);
          rn           = cf * rn;
          ain          = cf * ain;
          s1           = s1 + rn;
          s2           = s2 + ain;
          if ((fabs(rn) + fabs(ain)) < tol * (fabs(s1) + fabs(s2))) break;
          un  = un * ehn * ef;
          ehn = ehn * eh;
          en  = en * ey;
          vn  = vn + fn + fn + 1.0;
          fnh = fnh + 0.5;
          fn  = fn + 1.0;
        }
        s1 = s1 + s1;
        s2 = s2 + s2;
        if (z.real() == 0.0)
          s2 = s2 + yp;
        else {
          s1 = s1 + sxyh * sxyh / xp;
          s2 = s2 + sxy / tx;
        }
        // Power series for erf(xp), 0<=xp<=1
        RealType w  = 1.0;
        RealType ak = 1.5;
        RealType tm = 1.0;
        for (int i = 1; i <= 17; i++) {
          tm = tm * x2 / ak;
          w  = w + tm;
          if (tm <= tol) break;
          ak = ak + 1.0;
        }
        RealType ex = exp(-x2);
        w           = w * xp * fnorm * ex;
        RealType cf = ex / pi;
        s1          = cf * s1 + w;
        s2          = cf * s2;
        cans        = CmplxType(s1, s2);
        if (z.real() < 0.0) cans = -cans;
      }       // end (abs(yp) < 6.0)
      else {  //(abs(YP)>=6.0)
        // Asymptotic expansion for 0<=xp<=1 and abs(yp)>=6
        CmplxType rcz   = 0.5 / cz;
        CmplxType accum = CmplxType(1.0, 0.0);
        CmplxType term  = accum;
        RealType ak     = 1.0;
        for (int i = 1; i <= 35; i++) {
          term  = -term * ak * rcz;
          accum = accum + term;
          if (Kokkos::abs(term) / Kokkos::abs(accum) <= tol) break;
          ak = ak + 2.0;
        }
        accum       = accum * gnorm / zp;
        cz          = -cz;
        RealType er = cz.real();
        if (fabs(er) > 670.0) return CmplxType(inf, inf);
        RealType ei = cz.imag();
        cz          = exp(er) * CmplxType(cos(ei), sin(ei));
        cans        = 1.0 - accum * cz;
        if (z.real() < 0.0) cans = -cans;
      }  // end (abs(YP)>=6.0)
    }    // end (xp <= 1.0)
  }      // end (az > 2.0)
  return cans;
}

//! Compute scaled complementary error function erfcx(z)=exp(z^2)*erfc(z)
//! for z=cmplx(x,y).
template <class RealType>
KOKKOS_INLINE_FUNCTION Kokkos::complex<RealType> erfcx(
    const Kokkos::complex<RealType>& z) {
  // This function is a conversion of the corresponding Fortran program written
  // by D.E. Amos, May,1974. D.E. Amos' revisions of Jan 86 incorporated by
  // Ken Damrau on 27-Jan-1986 14:37:13
  //
  // Reference: NBS HANDBOOK OF MATHEMATICAL FUNCTIONS, AMS 55, By
  //           M. ABRAMOWITZ AND I.A. STEGUN, December,1955.
  // Summary:
  //  If x < 0, z is replaced by -z and all computation is done in the right
  //  half lane, except for z inside the circle abs(z)<=2, since
  //  erfc(-z)=2-erfc(z). The regions for computation are divided as follows
  //      (1)  abs(z)<=2 - Power series, NBS Handbook, p. 298
  //      (2)  abs(z)>2 and x>1 - continued fraction, NBS Handbook, p. 298
  //      (3)  abs(z)>2 and 0<=x<=1 and abs(y)<6 - series, NBS Handbook, p. 299
  //      (4)  abs(z)>2 and 0<=x<=1 and abs(y)>=6 - asymptotic expansion
  // Error condition: abs(z^2) > 670 is a fatal overflow error when x<0
  using Kokkos::cos;
  using Kokkos::exp;
  using Kokkos::fabs;
  using Kokkos::isinf;
  using Kokkos::sin;
  using Kokkos::Experimental::epsilon_v;
  using Kokkos::Experimental::infinity_v;
  using Kokkos::numbers::inv_sqrtpi_v;
  using Kokkos::numbers::pi_v;

  using CmplxType = Kokkos::complex<RealType>;

  constexpr auto inf = infinity_v<RealType>;
  constexpr auto tol = epsilon_v<RealType>;

  const RealType fnorm = 1.12837916709551;
  constexpr auto gnorm = inv_sqrtpi_v<RealType>;
  const RealType eh    = 0.606530659712633;
  const RealType ef    = 0.778800783071405;
  // const RealType tol   = 1.0e-13;
  constexpr auto pi = pi_v<RealType>;

  CmplxType cans;

  if ((isinf(z.real())) && (z.real() > 0)) {
    cans = CmplxType(0.0, 0.0);
    return cans;
  }
  if ((isinf(z.real())) && (z.real() < 0)) {
    cans = CmplxType(inf, inf);
    return cans;
  }

  RealType az = Kokkos::abs(z);
  if (az <= 2.0) {  // Series for abs(z)<=2.0
    CmplxType cz    = z * z;
    CmplxType accum = CmplxType(1.0, 0.0);
    CmplxType term  = accum;
    RealType ak     = 1.5;
    for (int i = 1; i <= 35; i++) {
      term  = term * cz / ak;
      accum = accum + term;
      if (Kokkos::abs(term) <= tol) break;
      ak = ak + 1.0;
    }
    cz          = -cz;
    RealType er = cz.real();
    RealType ei = cz.imag();
    accum       = accum * z * fnorm;
    cz          = exp(er) * CmplxType(cos(ei), sin(ei));
    cans        = 1.0 / cz - accum;
  }       // end (az <= 2.0)
  else {  //(az > 2.0)
    CmplxType zp = z;
    if (z.real() < 0.0) zp = -z;
    CmplxType cz = zp * zp;
    RealType xp  = zp.real();
    RealType yp  = zp.imag();
    if (xp > 1.0) {
      // continued fraction for erfc(z), abs(z)>2
      int n          = static_cast<int>(100.0 / az + 5.0);
      int fn         = n;
      CmplxType term = cz;
      for (int i = 1; i <= n; i++) {
        RealType fnh = fn - 0.5;
        term         = cz + (fnh * term) / (fn + term);
        fn           = fn - 1;
      }
      cans = zp * gnorm / term;
      if (z.real() >= 0.0) return cans;
      if (Kokkos::abs(cz) > 670.0) return CmplxType(inf, inf);
      ;
      cz          = -cz;
      RealType er = cz.real();
      RealType ei = cz.imag();
      cz          = exp(er) * CmplxType(cos(ei), sin(ei));
      cz          = 1.0 / cz;
      cans        = cz + cz - cans;
    }       // end (xp > 1.0)
    else {  //(xp <= 1.0)
      if (fabs(yp) <
          6.0) {  // Series (3) for abs(z)>2 and 0<=xp<=1 and abs(yp)<6
        RealType s1   = 0.0;
        RealType s2   = 0.0;
        RealType x2   = xp * xp;
        RealType fx2  = 4.0 * x2;
        RealType tx   = xp + xp;
        RealType xy   = xp * yp;
        RealType sxyh = sin(xy);
        RealType sxy  = sin(xy + xy);
        RealType cxy  = cos(xy + xy);
        RealType fn   = 1.0;
        RealType fnh  = 0.5;
        RealType ey   = exp(yp);
        RealType en   = ey;
        RealType ehn  = eh;
        RealType un   = ef;
        RealType vn   = 1.0;
        for (int i = 1; i <= 50; i++) {
          RealType ren = 1.0 / en;
          RealType csh = en + ren;
          RealType tm  = xp * csh;
          RealType ssh = en - ren;
          RealType tmp = fnh * ssh;
          RealType rn  = tx - tm * cxy + tmp * sxy;
          RealType ain = tm * sxy + tmp * cxy;
          RealType cf  = un / (vn + fx2);
          rn           = cf * rn;
          ain          = cf * ain;
          s1           = s1 + rn;
          s2           = s2 + ain;
          if ((fabs(rn) + fabs(ain)) < tol * (fabs(s1) + fabs(s2))) break;
          un  = un * ehn * ef;
          ehn = ehn * eh;
          en  = en * ey;
          vn  = vn + fn + fn + 1.0;
          fnh = fnh + 0.5;
          fn  = fn + 1.0;
        }
        s1 = s1 + s1;
        s2 = s2 + s2;
        if (z.real() == 0.0)
          s2 = s2 + yp;
        else {
          s1 = s1 + sxyh * sxyh / xp;
          s2 = s2 + sxy / tx;
        }
        // Power series for erf(xp), 0<=xp<=1
        RealType w  = 1.0;
        RealType ak = 1.5;
        RealType tm = 1.0;
        for (int i = 1; i <= 17; i++) {
          tm = tm * x2 / ak;
          w  = w + tm;
          if (tm <= tol) break;
          ak = ak + 1.0;
        }
        RealType ex   = exp(-x2);
        w             = w * xp * fnorm * ex;
        CmplxType rcz = CmplxType(cxy, sxy);
        RealType y2   = yp * yp;
        cz            = exp(x2 - y2) * rcz;
        rcz           = exp(-y2) * rcz;
        if (z.real() >= 0.0)
          cans = cz * (1.0 - w) - rcz * CmplxType(s1, s2) / pi;
        else
          cans = cz * (1.0 + w) + rcz * CmplxType(s1, s2) / pi;
      }       // end (abs(yp) < 6.0)
      else {  //(abs(YP)>=6.0)
        // Asymptotic expansion for 0<=xp<=1 and abs(yp)>=6
        CmplxType rcz   = 0.5 / cz;
        CmplxType accum = CmplxType(1.0, 0.0);
        CmplxType term  = accum;
        RealType ak     = 1.0;
        for (int i = 1; i <= 35; i++) {
          term  = -term * ak * rcz;
          accum = accum + term;
          if (Kokkos::abs(term) / Kokkos::abs(accum) <= tol) break;
          ak = ak + 2.0;
        }
        accum = accum * gnorm / zp;
        if (z.real() < 0.0) accum = -accum;
        cans = accum;
      }  // end (abs(YP)>=6.0)
    }    // end (xp <= 1.0)
  }      // end (az > 2.0)
  return cans;
}

//! Compute scaled complementary error function erfcx(x)=exp(x^2)*erfc(x)
//! for real x
template <class RealType>
KOKKOS_INLINE_FUNCTION RealType erfcx(RealType x) {
  using CmplxType = Kokkos::complex<RealType>;
  // Note: using erfcx(complex) for now
  // TODO: replace with an implementation of erfcx(real)
  CmplxType zin  = CmplxType(x, 0.0);
  CmplxType zout = erfcx(zin);
  return zout.real();
}

//! Compute Bessel function J0(z) of the first kind of order zero
//! for a complex argument
template <class CmplxType, class RealType, class IntType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_j0(const CmplxType& z,
                                               const RealType& joint_val = 25,
                                               const IntType& bw_start   = 70) {
  // This function is converted and modified from the corresponding Fortran
  // program CJYNB in S. Zhang & J. Jin "Computation of Special Functions"
  //(Wiley, 1996).
  // Input :  z         --- Complex argument
  //         joint_val --- Joint point of abs(z) separating small and large
  //                       argument regions
  //         bw_start  --- Starting point for backward recurrence
  // Output:  cbj0      --- J0(z)
  using Kokkos::fabs;
  using Kokkos::pow;
  using Kokkos::numbers::pi_v;

  CmplxType cbj0;
  constexpr auto pi    = pi_v<RealType>;
  const RealType a[12] = {
      -0.703125e-01,           0.112152099609375e+00,   -0.5725014209747314e+00,
      0.6074042001273483e+01,  -0.1100171402692467e+03, 0.3038090510922384e+04,
      -0.1188384262567832e+06, 0.6252951493434797e+07,  -0.4259392165047669e+09,
      0.3646840080706556e+11,  -0.3833534661393944e+13, 0.4854014686852901e+15};
  const RealType b[12] = {0.732421875e-01,        -0.2271080017089844e+00,
                          0.1727727502584457e+01, -0.2438052969955606e+02,
                          0.5513358961220206e+03, -0.1825775547429318e+05,
                          0.8328593040162893e+06, -0.5006958953198893e+08,
                          0.3836255180230433e+10, -0.3649010818849833e+12,
                          0.4218971570284096e+14, -0.5827244631566907e+16};

  RealType r2p = 2.0 / pi;
  RealType a0  = Kokkos::abs(z);
  RealType y0  = fabs(z.imag());
  CmplxType z1 = z;

  if (a0 < 1e-100) {  // Treat z=0 as a special case
    cbj0 = CmplxType(1.0, 0.0);
  } else {
    if (z.real() < 0.0) z1 = -z;
    if (a0 <= joint_val) {  // Using backward recurrence for |z|<=joint_val
                            // (default:25)
      CmplxType cbs = CmplxType(0.0, 0.0);
      CmplxType csu = CmplxType(0.0, 0.0);
      CmplxType csv = CmplxType(0.0, 0.0);
      CmplxType cf2 = CmplxType(0.0, 0.0);
      CmplxType cf1 = CmplxType(1e-100, 0.0);
      CmplxType cf, cs0;
      for (int k = bw_start; k >= 0; k--) {  // Backward recurrence (default:
                                             // 70)
        cf                    = 2.0 * (k + 1.0) / z * cf1 - cf2;
        RealType tmp_exponent = static_cast<RealType>(k / 2);
        if (k == 0) cbj0 = cf;
        if ((k == 2 * (k / 2)) && (k != 0)) {
          if (y0 <= 1.0)
            cbs = cbs + 2.0 * cf;
          else
            cbs = cbs + pow(-1.0, tmp_exponent) * 2.0 * cf;
          csu = csu + pow(-1.0, tmp_exponent) * cf / k;
        } else if (k > 1) {
          csv = csv + pow(-1.0, tmp_exponent) * k / (k * k - 1.0) * cf;
        }
        cf2 = cf1;
        cf1 = cf;
      }
      if (y0 <= 1.0)
        cs0 = cbs + cf;
      else
        cs0 = (cbs + cf) / Kokkos::cos(z);
      cbj0 = cbj0 / cs0;
    } else {  // Using asymptotic expansion (5.2.5) for |z|>joint_val
              // (default:25)
      CmplxType ct1 = z1 - 0.25 * pi;
      CmplxType cp0 = CmplxType(1.0, 0.0);
      for (int k = 1; k <= 12; k++) {  // Calculate (5.2.9)
        cp0 = cp0 + a[k - 1] * Kokkos::pow(z1, -2.0 * k);
      }
      CmplxType cq0 = -0.125 / z1;
      for (int k = 1; k <= 12; k++) {  // Calculate (5.2.10)
        cq0 = cq0 + b[k - 1] * Kokkos::pow(z1, -2.0 * k - 1);
      }
      CmplxType cu = Kokkos::sqrt(r2p / z1);
      cbj0         = cu * (cp0 * Kokkos::cos(ct1) - cq0 * Kokkos::sin(ct1));
    }
  }
  return cbj0;
}

//! Compute Bessel function Y0(z) of the second kind of order zero
//! for a complex argument
template <class CmplxType, class RealType, class IntType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_y0(const CmplxType& z,
                                               const RealType& joint_val = 25,
                                               const IntType& bw_start   = 70) {
  // This function is converted and modified from the corresponding Fortran
  // program CJYNB in S. Zhang & J. Jin "Computation of Special Functions"
  //(Wiley, 1996).
  //    Input :  z         --- Complex argument
  //             joint_val --- Joint point of abs(z) separating small and large
  //                           argument regions
  //             bw_start  --- Starting point for backward recurrence
  //    Output:  cby0      --- Y0(z)
  using Kokkos::fabs;
  using Kokkos::pow;
  using Kokkos::Experimental::infinity_v;
  using Kokkos::numbers::egamma_v;
  using Kokkos::numbers::pi_v;

  constexpr auto inf = infinity_v<RealType>;

  CmplxType cby0, cbj0;
  constexpr auto pi    = pi_v<RealType>;
  constexpr auto el    = egamma_v<RealType>;
  const RealType a[12] = {
      -0.703125e-01,           0.112152099609375e+00,   -0.5725014209747314e+00,
      0.6074042001273483e+01,  -0.1100171402692467e+03, 0.3038090510922384e+04,
      -0.1188384262567832e+06, 0.6252951493434797e+07,  -0.4259392165047669e+09,
      0.3646840080706556e+11,  -0.3833534661393944e+13, 0.4854014686852901e+15};
  const RealType b[12] = {0.732421875e-01,        -0.2271080017089844e+00,
                          0.1727727502584457e+01, -0.2438052969955606e+02,
                          0.5513358961220206e+03, -0.1825775547429318e+05,
                          0.8328593040162893e+06, -0.5006958953198893e+08,
                          0.3836255180230433e+10, -0.3649010818849833e+12,
                          0.4218971570284096e+14, -0.5827244631566907e+16};

  RealType r2p = 2.0 / pi;
  RealType a0  = Kokkos::abs(z);
  RealType y0  = fabs(z.imag());
  CmplxType ci = CmplxType(0.0, 1.0);
  CmplxType z1 = z;

  if (a0 < 1e-100) {  // Treat z=0 as a special case
    cby0 = -CmplxType(inf, 0.0);
  } else {
    if (z.real() < 0.0) z1 = -z;
    if (a0 <= joint_val) {  // Using backward recurrence for |z|<=joint_val
                            // (default:25)
      CmplxType cbs = CmplxType(0.0, 0.0);
      CmplxType csu = CmplxType(0.0, 0.0);
      CmplxType csv = CmplxType(0.0, 0.0);
      CmplxType cf2 = CmplxType(0.0, 0.0);
      CmplxType cf1 = CmplxType(1e-100, 0.0);
      CmplxType cf, cs0, ce;
      for (int k = bw_start; k >= 0; k--) {  // Backward recurrence (default:
                                             // 70)
        cf                    = 2.0 * (k + 1.0) / z * cf1 - cf2;
        RealType tmp_exponent = static_cast<RealType>(k / 2);
        if (k == 0) cbj0 = cf;
        if ((k == 2 * (k / 2)) && (k != 0)) {
          if (y0 <= 1.0)
            cbs = cbs + 2.0 * cf;
          else
            cbs = cbs + pow(-1.0, tmp_exponent) * 2.0 * cf;
          csu = csu + pow(-1.0, tmp_exponent) * cf / k;
        } else if (k > 1) {
          csv = csv + pow(-1.0, tmp_exponent) * k / (k * k - 1.0) * cf;
        }
        cf2 = cf1;
        cf1 = cf;
      }
      if (y0 <= 1.0)
        cs0 = cbs + cf;
      else
        cs0 = (cbs + cf) / Kokkos::cos(z);
      cbj0 = cbj0 / cs0;
      ce   = Kokkos::log(z / 2.0) + el;
      cby0 = r2p * (ce * cbj0 - 4.0 * csu / cs0);
    } else {  // Using asymptotic expansion (5.2.6) for |z|>joint_val
              // (default:25)
      CmplxType ct1 = z1 - 0.25 * pi;
      CmplxType cp0 = CmplxType(1.0, 0.0);
      for (int k = 1; k <= 12; k++) {  // Calculate (5.2.9)
        cp0 = cp0 + a[k - 1] * Kokkos::pow(z1, -2.0 * k);
      }
      CmplxType cq0 = -0.125 / z1;
      for (int k = 1; k <= 12; k++) {  // Calculate (5.2.10)
        cq0 = cq0 + b[k - 1] * Kokkos::pow(z1, -2.0 * k - 1);
      }
      CmplxType cu = Kokkos::sqrt(r2p / z1);
      cbj0         = cu * (cp0 * Kokkos::cos(ct1) - cq0 * Kokkos::sin(ct1));
      cby0         = cu * (cp0 * Kokkos::sin(ct1) + cq0 * Kokkos::cos(ct1));

      if (z.real() < 0.0) {  // Apply (5.4.2)
        if (z.imag() < 0.0) cby0 = cby0 - 2.0 * ci * cbj0;
        if (z.imag() >= 0.0) cby0 = cby0 + 2.0 * ci * cbj0;
      }
    }
  }
  return cby0;
}

//! Compute Bessel function J1(z) of the first kind of order one
//! for a complex argument
template <class CmplxType, class RealType, class IntType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_j1(const CmplxType& z,
                                               const RealType& joint_val = 25,
                                               const IntType& bw_start   = 70) {
  // This function is converted and modified from the corresponding Fortran
  // program CJYNB in S. Zhang & J. Jin "Computation of Special Functions"
  //(Wiley, 1996).
  //    Input :  z         --- Complex argument
  //             joint_val --- Joint point of abs(z) separating small and large
  //                           argument regions
  //             bw_start  --- Starting point for backward recurrence
  //    Output:  cbj1      --- J1(z)
  using Kokkos::fabs;
  using Kokkos::pow;
  using Kokkos::numbers::pi_v;

  CmplxType cbj1;
  constexpr auto pi     = pi_v<RealType>;
  const RealType a1[12] = {0.1171875e+00,          -0.144195556640625e+00,
                           0.6765925884246826e+00, -0.6883914268109947e+01,
                           0.1215978918765359e+03, -0.3302272294480852e+04,
                           0.1276412726461746e+06, -0.6656367718817688e+07,
                           0.4502786003050393e+09, -0.3833857520742790e+11,
                           0.4011838599133198e+13, -0.5060568503314727e+15};
  const RealType b1[12] = {
      -0.1025390625e+00,       0.2775764465332031e+00,  -0.1993531733751297e+01,
      0.2724882731126854e+02,  -0.6038440767050702e+03, 0.1971837591223663e+05,
      -0.8902978767070678e+06, 0.5310411010968522e+08,  -0.4043620325107754e+10,
      0.3827011346598605e+12,  -0.4406481417852278e+14, 0.6065091351222699e+16};

  RealType r2p = 2.0 / pi;
  RealType a0  = Kokkos::abs(z);
  RealType y0  = fabs(z.imag());
  CmplxType z1 = z;

  if (a0 < 1e-100) {  // Treat z=0 as a special case
    cbj1 = CmplxType(0.0, 0.0);
  } else {
    if (z.real() < 0.0) z1 = -z;
    if (a0 <= joint_val) {  // Using backward recurrence for |z|<=joint_val
                            // (default:25)
      CmplxType cbs = CmplxType(0.0, 0.0);
      CmplxType csu = CmplxType(0.0, 0.0);
      CmplxType csv = CmplxType(0.0, 0.0);
      CmplxType cf2 = CmplxType(0.0, 0.0);
      CmplxType cf1 = CmplxType(1e-100, 0.0);
      CmplxType cf, cs0;
      for (int k = bw_start; k >= 0; k--) {  // Backward recurrence (default:
                                             // 70)
        cf                    = 2.0 * (k + 1.0) / z * cf1 - cf2;
        RealType tmp_exponent = static_cast<RealType>(k / 2);
        if (k == 1) cbj1 = cf;
        if ((k == 2 * (k / 2)) && (k != 0)) {
          if (y0 <= 1.0)
            cbs = cbs + 2.0 * cf;
          else
            cbs = cbs + pow(-1.0, tmp_exponent) * 2.0 * cf;
          csu = csu + pow(-1.0, tmp_exponent) * cf / k;
        } else if (k > 1) {
          csv = csv + pow(-1.0, tmp_exponent) * k / (k * k - 1.0) * cf;
        }
        cf2 = cf1;
        cf1 = cf;
      }
      if (y0 <= 1.0)
        cs0 = cbs + cf;
      else
        cs0 = (cbs + cf) / Kokkos::cos(z);
      cbj1 = cbj1 / cs0;
    } else {  // Using asymptotic expansion (5.2.5) for |z|>joint_val
              // (default:25)
      CmplxType ct2 = z1 - 0.75 * pi;
      CmplxType cp1 = CmplxType(1.0, 0.0);
      for (int k = 1; k <= 12; k++) {  // Calculate (5.2.11)
        cp1 = cp1 + a1[k - 1] * Kokkos::pow(z1, -2.0 * k);
      }
      CmplxType cq1 = 0.375 / z1;
      for (int k = 1; k <= 12; k++) {  // Calculate (5.2.12)
        cq1 = cq1 + b1[k - 1] * Kokkos::pow(z1, -2.0 * k - 1);
      }
      CmplxType cu = Kokkos::sqrt(r2p / z1);
      cbj1         = cu * (cp1 * Kokkos::cos(ct2) - cq1 * Kokkos::sin(ct2));

      if (real(z) < 0.0) {  // Apply (5.4.2)
        cbj1 = -cbj1;
      }
    }
  }
  return cbj1;
}

//! Compute Bessel function Y1(z) of the second kind of order one
//! for a complex argument
template <class CmplxType, class RealType, class IntType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_y1(const CmplxType& z,
                                               const RealType& joint_val = 25,
                                               const IntType& bw_start   = 70) {
  // This function is converted and modified from the corresponding Fortran
  // program CJYNB in S. Zhang & J. Jin "Computation of Special Functions"
  //(Wiley, 1996).
  //    Input :  z         --- Complex argument
  //             joint_val --- Joint point of abs(z) separating small and large
  //                           argument regions
  //             bw_start  --- Starting point for backward recurrence
  //    Output:  cby1      --- Y1(z)
  using Kokkos::fabs;
  using Kokkos::pow;
  using Kokkos::Experimental::infinity_v;
  using Kokkos::numbers::egamma_v;
  using Kokkos::numbers::pi_v;

  constexpr auto inf = infinity_v<RealType>;

  CmplxType cby1, cbj0, cbj1, cby0;
  constexpr auto pi     = pi_v<RealType>;
  constexpr auto el     = egamma_v<RealType>;
  const RealType a1[12] = {0.1171875e+00,          -0.144195556640625e+00,
                           0.6765925884246826e+00, -0.6883914268109947e+01,
                           0.1215978918765359e+03, -0.3302272294480852e+04,
                           0.1276412726461746e+06, -0.6656367718817688e+07,
                           0.4502786003050393e+09, -0.3833857520742790e+11,
                           0.4011838599133198e+13, -0.5060568503314727e+15};
  const RealType b1[12] = {
      -0.1025390625e+00,       0.2775764465332031e+00,  -0.1993531733751297e+01,
      0.2724882731126854e+02,  -0.6038440767050702e+03, 0.1971837591223663e+05,
      -0.8902978767070678e+06, 0.5310411010968522e+08,  -0.4043620325107754e+10,
      0.3827011346598605e+12,  -0.4406481417852278e+14, 0.6065091351222699e+16};

  RealType r2p = 2.0 / pi;
  RealType a0  = Kokkos::abs(z);
  RealType y0  = fabs(z.imag());
  CmplxType ci = CmplxType(0.0, 1.0);
  CmplxType z1 = z;

  if (a0 < 1e-100) {  // Treat z=0 as a special case
    cby1 = -CmplxType(inf, 0.0);
  } else {
    if (z.real() < 0.0) z1 = -z;
    if (a0 <= joint_val) {  // Using backward recurrence for |z|<=joint_val
                            // (default:25)
      CmplxType cbs = CmplxType(0.0, 0.0);
      CmplxType csu = CmplxType(0.0, 0.0);
      CmplxType csv = CmplxType(0.0, 0.0);
      CmplxType cf2 = CmplxType(0.0, 0.0);
      CmplxType cf1 = CmplxType(1e-100, 0.0);
      CmplxType cf, cs0, ce;
      for (int k = bw_start; k >= 0; k--) {  // Backward recurrence (default:
                                             // 70)
        cf                    = 2.0 * (k + 1.0) / z * cf1 - cf2;
        RealType tmp_exponent = static_cast<RealType>(k / 2);
        if (k == 1) cbj1 = cf;
        if (k == 0) cbj0 = cf;
        if ((k == 2 * (k / 2)) && (k != 0)) {
          if (y0 <= 1.0)
            cbs = cbs + 2.0 * cf;
          else
            cbs = cbs + pow(-1.0, tmp_exponent) * 2.0 * cf;
          csu = csu + pow(-1.0, tmp_exponent) * cf / k;
        } else if (k > 1) {
          csv = csv + pow(-1.0, tmp_exponent) * k / (k * k - 1.0) * cf;
        }
        cf2 = cf1;
        cf1 = cf;
      }
      if (y0 <= 1.0)
        cs0 = cbs + cf;
      else
        cs0 = (cbs + cf) / Kokkos::cos(z);
      cbj0 = cbj0 / cs0;
      ce   = Kokkos::log(z / 2.0) + el;
      cby0 = r2p * (ce * cbj0 - 4.0 * csu / cs0);
      cbj1 = cbj1 / cs0;
      cby1 = (cbj1 * cby0 - 2.0 / (pi * z)) / cbj0;
    } else {  // Using asymptotic expansion (5.2.5) for |z|>joint_val
              // (default:25)
      CmplxType ct2 = z1 - 0.75 * pi;
      CmplxType cp1 = CmplxType(1.0, 0.0);
      for (int k = 1; k <= 12; k++) {  // Calculate (5.2.11)
        cp1 = cp1 + a1[k - 1] * Kokkos::pow(z1, -2.0 * k);
      }
      CmplxType cq1 = 0.375 / z1;
      for (int k = 1; k <= 12; k++) {  // Calculate (5.2.12)
        cq1 = cq1 + b1[k - 1] * Kokkos::pow(z1, -2.0 * k - 1);
      }
      CmplxType cu = Kokkos::sqrt(r2p / z1);
      cbj1         = cu * (cp1 * Kokkos::cos(ct2) - cq1 * Kokkos::sin(ct2));
      cby1         = cu * (cp1 * Kokkos::sin(ct2) + cq1 * Kokkos::cos(ct2));

      if (z.real() < 0.0) {  // Apply (5.4.2)
        if (z.imag() < 0.0) cby1 = -(cby1 - 2.0 * ci * cbj1);
        if (z.imag() >= 0.0) cby1 = -(cby1 + 2.0 * ci * cbj1);
      }
    }
  }
  return cby1;
}

//! Compute modified Bessel function I0(z) of the first kind of order zero
//! for a complex argument
template <class CmplxType, class RealType, class IntType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_i0(const CmplxType& z,
                                               const RealType& joint_val = 25,
                                               const IntType& bw_start   = 70) {
  // This function is converted and modified from the corresponding Fortran
  // programs CIKNB and CIK01 in S. Zhang & J. Jin "Computation of Special
  // Functions" (Wiley, 1996).
  //    Input :  z         --- Complex argument
  //             joint_val --- Joint point of abs(z) separating small and large
  //                           argument regions
  //             bw_start  --- Starting point for backward recurrence
  //    Output:  cbi0      --- I0(z)
  using Kokkos::numbers::pi_v;

  CmplxType cbi0;
  constexpr auto pi    = pi_v<RealType>;
  const RealType a[12] = {0.125,
                          7.03125e-2,
                          7.32421875e-2,
                          1.1215209960938e-1,
                          2.2710800170898e-1,
                          5.7250142097473e-1,
                          1.7277275025845e0,
                          6.0740420012735e0,
                          2.4380529699556e1,
                          1.1001714026925e2,
                          5.5133589612202e2,
                          3.0380905109224e3};

  RealType a0  = Kokkos::abs(z);
  CmplxType z1 = z;

  if (a0 < 1e-100) {  // Treat z=0 as a special case
    cbi0 = CmplxType(1.0, 0.0);
  } else {
    if (z.real() < 0.0) z1 = -z;
    if (a0 <= joint_val) {  // Using backward recurrence for |z|<=joint_val
                            // (default:25)
      CmplxType cbs = CmplxType(0.0, 0.0);
      // CmplxType csk0 = CmplxType(0.0,0.0);
      CmplxType cf0 = CmplxType(0.0, 0.0);
      CmplxType cf1 = CmplxType(1e-100, 0.0);
      CmplxType cf, cs0;
      for (int k = bw_start; k >= 0; k--) {  // Backward recurrence (default:
                                             // 70)
        cf = 2.0 * (k + 1.0) * cf1 / z1 + cf0;
        if (k == 0) cbi0 = cf;
        // if ((k == 2*(k/2)) && (k != 0)) {
        //  csk0 = csk0+4.0*cf/static_cast<RealType>(k);
        //}
        cbs = cbs + 2.0 * cf;
        cf0 = cf1;
        cf1 = cf;
      }
      cs0  = Kokkos::exp(z1) / (cbs - cf);
      cbi0 = cbi0 * cs0;
    } else {  // Using asymptotic expansion (6.2.1) for |z|>joint_val
              // (default:25)
      CmplxType ca = Kokkos::exp(z1) / Kokkos::sqrt(2.0 * pi * z1);
      cbi0         = CmplxType(1.0, 0.0);
      CmplxType zr = 1.0 / z1;
      for (int k = 1; k <= 12; k++) {
        cbi0 = cbi0 + a[k - 1] * Kokkos::pow(zr, 1.0 * k);
      }
      cbi0 = ca * cbi0;
    }
  }
  return cbi0;
}

//! Compute modified Bessel function K0(z) of the second kind of order zero
//! for a complex argument
template <class CmplxType, class RealType, class IntType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_k0(const CmplxType& z,
                                               const RealType& joint_val = 9,
                                               const IntType& bw_start   = 30) {
  // This function is converted and modified from the corresponding Fortran
  // programs CIKNB and CIK01 in S. Zhang & J. Jin "Computation of Special
  // Functions" (Wiley, 1996).
  //    Purpose: Compute modified Bessel function K0(z) of the second kind of
  //             order zero for a complex argument
  //    Input :  z         --- Complex argument
  //             joint_val --- Joint point of abs(z) separating small and large
  //                           argument regions
  //             bw_start  --- Starting point for backward recurrence
  //    Output:  cbk0      --- K0(z)
  using Kokkos::pow;
  using Kokkos::Experimental::infinity_v;
  using Kokkos::numbers::egamma_v;
  using Kokkos::numbers::pi_v;

  constexpr auto inf = infinity_v<RealType>;

  CmplxType cbk0, cbi0;
  constexpr auto pi = pi_v<RealType>;
  constexpr auto el = egamma_v<RealType>;

  RealType a0  = Kokkos::abs(z);
  CmplxType ci = CmplxType(0.0, 1.0);
  CmplxType z1 = z;

  if (a0 < 1e-100) {  // Treat z=0 as a special case
    cbk0 = CmplxType(inf, 0.0);
  } else {
    if (z.real() < 0.0) z1 = -z;
    if (a0 <= joint_val) {  // Using backward recurrence for |z|<=joint_val
                            // (default:9)
      CmplxType cbs  = CmplxType(0.0, 0.0);
      CmplxType csk0 = CmplxType(0.0, 0.0);
      CmplxType cf0  = CmplxType(0.0, 0.0);
      CmplxType cf1  = CmplxType(1e-100, 0.0);
      CmplxType cf, cs0;
      for (int k = bw_start; k >= 0; k--) {  // Backward recurrence (default:
                                             // 30)
        cf = 2.0 * (k + 1.0) * cf1 / z1 + cf0;
        if (k == 0) cbi0 = cf;
        if ((k == 2 * (k / 2)) && (k != 0)) {
          csk0 = csk0 + 4.0 * cf / static_cast<RealType>(k);
        }
        cbs = cbs + 2.0 * cf;
        cf0 = cf1;
        cf1 = cf;
      }
      cs0  = Kokkos::exp(z1) / (cbs - cf);
      cbi0 = cbi0 * cs0;
      cbk0 = -(Kokkos::log(0.5 * z1) + el) * cbi0 + cs0 * csk0;
    } else {  // Using asymptotic expansion (6.2.2) for |z|>joint_val
              // (default:9)
      CmplxType ca0  = Kokkos::sqrt(pi / (2.0 * z1)) * Kokkos::exp(-z1);
      CmplxType cbkl = CmplxType(1.0, 0.0);
      CmplxType cr   = CmplxType(1.0, 0.0);
      for (int k = 1; k <= 30; k++) {
        cr   = 0.125 * cr * (0.0 - pow(2.0 * k - 1.0, 2.0)) / (k * z1);
        cbkl = cbkl + cr;
      }
      cbk0 = ca0 * cbkl;
    }
    if (z.real() < 0.0) {  // Apply (6.4.4)
      if (z.imag() < 0.0)
        cbk0 = cbk0 + ci * pi * cyl_bessel_i0<CmplxType, RealType, IntType>(z);
      if (z.imag() >= 0.0)
        cbk0 = cbk0 - ci * pi * cyl_bessel_i0<CmplxType, RealType, IntType>(z);
    }
  }
  return cbk0;
}

//! Compute modified Bessel function I1(z) of the first kind of order one
//! for a complex argument
template <class CmplxType, class RealType, class IntType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_i1(const CmplxType& z,
                                               const RealType& joint_val = 25,
                                               const IntType& bw_start   = 70) {
  // This function is converted and modified from the corresponding Fortran
  // programs CIKNB and CIK01 in S. Zhang & J. Jin "Computation of Special
  // Functions" (Wiley, 1996).
  //    Input :  z         --- Complex argument
  //             joint_val --- Joint point of abs(z) separating small and large
  //                           argument regions
  //             bw_start  --- Starting point for backward recurrence
  //    Output:  cbi1      --- I1(z)
  using Kokkos::numbers::pi_v;

  CmplxType cbi1;
  constexpr auto pi    = pi_v<RealType>;
  const RealType b[12] = {-0.375,
                          -1.171875e-1,
                          -1.025390625e-1,
                          -1.4419555664063e-1,
                          -2.7757644653320e-1,
                          -6.7659258842468e-1,
                          -1.9935317337513,
                          -6.8839142681099,
                          -2.7248827311269e1,
                          -1.2159789187654e2,
                          -6.0384407670507e2,
                          -3.3022722944809e3};

  RealType a0  = Kokkos::abs(z);
  CmplxType z1 = z;

  if (a0 < 1e-100) {  // Treat z=0 as a special case
    cbi1 = CmplxType(0.0, 0.0);
  } else {
    if (z.real() < 0.0) z1 = -z;
    if (a0 <= joint_val) {  // Using backward recurrence for |z|<=joint_val
                            // (default:25)
      CmplxType cbs = CmplxType(0.0, 0.0);
      // CmplxType csk0 = CmplxType(0.0,0.0);
      CmplxType cf0 = CmplxType(0.0, 0.0);
      CmplxType cf1 = CmplxType(1e-100, 0.0);
      CmplxType cf, cs0;
      for (int k = bw_start; k >= 0; k--) {  // Backward recurrence (default:
                                             // 70)
        cf = 2.0 * (k + 1.0) * cf1 / z1 + cf0;
        if (k == 1) cbi1 = cf;
        // if ((k == 2*(k/2)) && (k != 0)) {
        //  csk0 = csk0+4.0*cf/static_cast<RealType>(k);
        //}
        cbs = cbs + 2.0 * cf;
        cf0 = cf1;
        cf1 = cf;
      }
      cs0  = Kokkos::exp(z1) / (cbs - cf);
      cbi1 = cbi1 * cs0;
    } else {  // Using asymptotic expansion (6.2.1) for |z|>joint_val
              // (default:25)
      CmplxType ca = Kokkos::exp(z1) / Kokkos::sqrt(2.0 * pi * z1);
      cbi1         = CmplxType(1.0, 0.0);
      CmplxType zr = 1.0 / z1;
      for (int k = 1; k <= 12; k++) {
        cbi1 = cbi1 + b[k - 1] * Kokkos::pow(zr, 1.0 * k);
      }
      cbi1 = ca * cbi1;
    }
    if (z.real() < 0.0) {  // Apply (6.4.4)
      cbi1 = -cbi1;
    }
  }
  return cbi1;
}

//! Compute modified Bessel function K1(z) of the second kind of order one
//! for a complex argument
template <class CmplxType, class RealType, class IntType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_k1(const CmplxType& z,
                                               const RealType& joint_val = 9,
                                               const IntType& bw_start   = 30) {
  // This function is converted and modified from the corresponding Fortran
  // programs CIKNB and CIK01 in S. Zhang & J. Jin "Computation of Special
  // Functions" (Wiley, 1996).
  //    Input :  z         --- Complex argument
  //             joint_val --- Joint point of abs(z) separating small and large
  //                           argument regions
  //             bw_start  --- Starting point for backward recurrence
  //    Output:  cbk1      --- K1(z)
  using Kokkos::pow;
  using Kokkos::Experimental::infinity_v;
  using Kokkos::numbers::egamma_v;
  using Kokkos::numbers::pi_v;

  constexpr auto inf = infinity_v<RealType>;

  CmplxType cbk0, cbi0, cbk1, cbi1;
  constexpr auto pi = pi_v<RealType>;
  constexpr auto el = egamma_v<RealType>;

  RealType a0  = Kokkos::abs(z);
  CmplxType ci = CmplxType(0.0, 1.0);
  CmplxType z1 = z;

  if (a0 < 1e-100) {  // Treat z=0 as a special case
    cbk1 = CmplxType(inf, 0.0);
  } else {
    if (z.real() < 0.0) z1 = -z;
    if (a0 <= joint_val) {  // Using backward recurrence for |z|<=joint_val
                            // (default:9)
      CmplxType cbs  = CmplxType(0.0, 0.0);
      CmplxType csk0 = CmplxType(0.0, 0.0);
      CmplxType cf0  = CmplxType(0.0, 0.0);
      CmplxType cf1  = CmplxType(1e-100, 0.0);
      CmplxType cf, cs0;
      for (int k = bw_start; k >= 0; k--) {  // Backward recurrence (default:
                                             // 30)
        cf = 2.0 * (k + 1.0) * cf1 / z1 + cf0;
        if (k == 1) cbi1 = cf;
        if (k == 0) cbi0 = cf;
        if ((k == 2 * (k / 2)) && (k != 0)) {
          csk0 = csk0 + 4.0 * cf / static_cast<RealType>(k);
        }
        cbs = cbs + 2.0 * cf;
        cf0 = cf1;
        cf1 = cf;
      }
      cs0  = Kokkos::exp(z1) / (cbs - cf);
      cbi0 = cbi0 * cs0;
      cbi1 = cbi1 * cs0;
      cbk0 = -(Kokkos::log(0.5 * z1) + el) * cbi0 + cs0 * csk0;
      cbk1 = (1.0 / z1 - cbi1 * cbk0) / cbi0;
    } else {  // Using asymptotic expansion (6.2.2) for |z|>joint_val
              // (default:9)
      CmplxType ca0  = Kokkos::sqrt(pi / (2.0 * z1)) * Kokkos::exp(-z1);
      CmplxType cbkl = CmplxType(1.0, 0.0);
      CmplxType cr   = CmplxType(1.0, 0.0);
      for (int k = 1; k <= 30; k++) {
        cr   = 0.125 * cr * (4.0 - pow(2.0 * k - 1.0, 2.0)) / (k * z1);
        cbkl = cbkl + cr;
      }
      cbk1 = ca0 * cbkl;
    }
    if (z.real() < 0.0) {  // Apply (6.4.4)
      if (z.imag() < 0.0)
        cbk1 = -cbk1 - ci * pi * cyl_bessel_i1<CmplxType, RealType, IntType>(z);
      if (z.imag() >= 0.0)
        cbk1 = -cbk1 + ci * pi * cyl_bessel_i1<CmplxType, RealType, IntType>(z);
    }
  }
  return cbk1;
}

//! Compute Hankel function H10(z) of the first kind of order zero
//! for a complex argument
template <class CmplxType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_h10(const CmplxType& z) {
  // This function is converted and modified from the corresponding Fortran
  // programs CH12N in S. Zhang & J. Jin "Computation of Special Functions"
  //(Wiley, 1996).
  using RealType = typename CmplxType::value_type;
  using Kokkos::Experimental::infinity_v;
  using Kokkos::numbers::pi_v;

  constexpr auto inf = infinity_v<RealType>;

  CmplxType ch10, cbk0, cbj0, cby0;
  constexpr auto pi = pi_v<RealType>;
  CmplxType ci      = CmplxType(0.0, 1.0);

  if ((z.real() == 0.0) && (z.imag() == 0.0)) {
    ch10 = CmplxType(1.0, -inf);
  } else if (z.imag() <= 0.0) {
    cbj0 = cyl_bessel_j0<CmplxType, RealType, int>(z);
    cby0 = cyl_bessel_y0<CmplxType, RealType, int>(z);
    ch10 = cbj0 + ci * cby0;
  } else {  //(z.imag() > 0.0)
    cbk0 = cyl_bessel_k0<CmplxType, RealType, int>(-ci * z, 18.0, 70);
    ch10 = 2.0 / (pi * ci) * cbk0;
  }

  return ch10;
}

//! Compute Hankel function H11(z) of the first kind of order one
//! for a complex argument
template <class CmplxType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_h11(const CmplxType& z) {
  // This function is converted and modified from the corresponding Fortran
  // programs CH12N in S. Zhang & J. Jin "Computation of Special Functions"
  //(Wiley, 1996).
  using RealType = typename CmplxType::value_type;
  using Kokkos::Experimental::infinity_v;
  using Kokkos::numbers::pi_v;

  constexpr auto inf = infinity_v<RealType>;

  CmplxType ch11, cbk1, cbj1, cby1;
  constexpr auto pi = pi_v<RealType>;
  CmplxType ci      = CmplxType(0.0, 1.0);

  if ((z.real() == 0.0) && (z.imag() == 0.0)) {
    ch11 = CmplxType(0.0, -inf);
  } else if (z.imag() <= 0.0) {
    cbj1 = cyl_bessel_j1<CmplxType, RealType, int>(z);
    cby1 = cyl_bessel_y1<CmplxType, RealType, int>(z);
    ch11 = cbj1 + ci * cby1;
  } else {  //(z.imag() > 0.0)
    cbk1 = cyl_bessel_k1<CmplxType, RealType, int>(-ci * z, 18.0, 70);
    ch11 = -2.0 / pi * cbk1;
  }

  return ch11;
}

//! Compute Hankel function H20(z) of the second kind of order zero
//! for a complex argument
template <class CmplxType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_h20(const CmplxType& z) {
  // This function is converted and modified from the corresponding Fortran
  // programs CH12N in S. Zhang & J. Jin "Computation of Special Functions"
  //(Wiley, 1996).
  using RealType = typename CmplxType::value_type;
  using Kokkos::Experimental::infinity_v;
  using Kokkos::numbers::pi_v;

  constexpr auto inf = infinity_v<RealType>;

  CmplxType ch20, cbk0, cbj0, cby0;
  constexpr auto pi = pi_v<RealType>;
  CmplxType ci      = CmplxType(0.0, 1.0);

  if ((z.real() == 0.0) && (z.imag() == 0.0)) {
    ch20 = CmplxType(1.0, inf);
  } else if (z.imag() >= 0.0) {
    cbj0 = cyl_bessel_j0<CmplxType, RealType, int>(z);
    cby0 = cyl_bessel_y0<CmplxType, RealType, int>(z);
    ch20 = cbj0 - ci * cby0;
  } else {  //(z.imag() < 0.0)
    cbk0 = cyl_bessel_k0<CmplxType, RealType, int>(ci * z, 18.0, 70);
    ch20 = 2.0 / pi * ci * cbk0;
  }

  return ch20;
}

//! Compute Hankel function H21(z) of the second kind of order one
//! for a complex argument
template <class CmplxType>
KOKKOS_INLINE_FUNCTION CmplxType cyl_bessel_h21(const CmplxType& z) {
  // This function is converted and modified from the corresponding Fortran
  // programs CH12N in S. Zhang & J. Jin "Computation of Special Functions"
  //(Wiley, 1996).
  using RealType = typename CmplxType::value_type;
  using Kokkos::Experimental::infinity_v;
  using Kokkos::numbers::pi_v;

  constexpr auto inf = infinity_v<RealType>;

  CmplxType ch21, cbk1, cbj1, cby1;
  constexpr auto pi = pi_v<RealType>;
  CmplxType ci      = CmplxType(0.0, 1.0);

  if ((z.real() == 0.0) && (z.imag() == 0.0)) {
    ch21 = CmplxType(0.0, inf);
  } else if (z.imag() >= 0.0) {
    cbj1 = cyl_bessel_j1<CmplxType, RealType, int>(z);
    cby1 = cyl_bessel_y1<CmplxType, RealType, int>(z);
    ch21 = cbj1 - ci * cby1;
  } else {  //(z.imag() < 0.0)
    cbk1 = cyl_bessel_k1<CmplxType, RealType, int>(ci * z, 18.0, 70);
    ch21 = -2.0 / pi * cbk1;
  }

  return ch21;
}

}  // namespace Experimental
}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_MATHSPECFUNCTIONS
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_MATHSPECFUNCTIONS
#endif
#endif

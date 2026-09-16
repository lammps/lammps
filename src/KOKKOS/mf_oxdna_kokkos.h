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

#ifndef MF_OXDNA_KOKKOS_H
#define MF_OXDNA_KOKKOS_H

#include "kokkos_type.h"

namespace MFOxdnaKokkos {

#if defined(KOKKOS_ENABLE_HIP)
template<class DeviceType, class Tag>
using OxdnaRangePolicy = Kokkos::RangePolicy<DeviceType, Tag, Kokkos::LaunchBounds<128, 1>>;
#elif defined(KOKKOS_ENABLE_CUDA)
template<class DeviceType, class Tag>
using OxdnaRangePolicy = Kokkos::RangePolicy<DeviceType, Tag, Kokkos::LaunchBounds<64, 1>>;
#else
template<class DeviceType, class Tag>
using OxdnaRangePolicy = Kokkos::RangePolicy<DeviceType, Tag>;
#endif

/* ----------------------------------------------------------------------
   f1 modulation factor
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT F1_KK(KK_FLOAT r, KK_FLOAT eps, KK_FLOAT a, KK_FLOAT cut_0,
                     KK_FLOAT cut_lc, KK_FLOAT cut_hc, KK_FLOAT cut_lo,
                     KK_FLOAT cut_hi, KK_FLOAT b_lo,
                     KK_FLOAT b_hi, KK_FLOAT shift)
{
  if (r > cut_hc) {
    return 0.0;
  } else if (r > cut_hi) {
    return eps * b_hi * (r - cut_hc) * (r - cut_hc);
  } else if (r > cut_lo) {
    KK_FLOAT tmp = 1 - Kokkos::expf(-(r - cut_0) * a);
    return eps * tmp * tmp - shift;
  } else if (r > cut_lc) {
    return eps * b_lo * (r - cut_lc) * (r - cut_lc);
  } else {
    return 0.0;
  }
}

KOKKOS_INLINE_FUNCTION
static KK_FLOAT F1_KK(KK_FLOAT r, KK_FLOAT eps, KK_FLOAT a, KK_FLOAT cut_0,
                     KK_FLOAT cut_lc, KK_FLOAT cut_hc, KK_FLOAT cut_lo,
                     KK_FLOAT cut_hi, KK_FLOAT b_lo,
                     KK_FLOAT b_hi, KK_FLOAT shift, KK_FLOAT &df1)
{
  if (r > cut_hc) {
    df1 = 0.0;
    return 0.0;
  } else if (r > cut_hi) {
    df1 = 2 * eps * b_hi * (1 - cut_hc / r);
    return eps * b_hi * (r - cut_hc) * (r - cut_hc);
  } else if (r > cut_lo) {
    KK_FLOAT tmp = Kokkos::expf(-(r - cut_0) * a);
    df1 = 2 * eps * (1 - tmp) * tmp * a / r;
    tmp = 1 - tmp;
    return eps * tmp * tmp - shift;
  } else if (r > cut_lc) {
    df1 = 2 * eps * b_lo * (1 - cut_lc / r);
    return eps * b_lo * (r - cut_lc) * (r - cut_lc);
  } else {
    df1 = 0.0;
    return 0.0;
  }
}

/* ----------------------------------------------------------------------
   derivative of f1 modulation factor
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT DF1_KK(KK_FLOAT r, KK_FLOAT eps, KK_FLOAT a, KK_FLOAT cut_0,
                      KK_FLOAT cut_lc, KK_FLOAT cut_hc, KK_FLOAT cut_lo,
                      KK_FLOAT cut_hi, KK_FLOAT b_lo, KK_FLOAT b_hi)
{
  if (r > cut_hc) {
    return 0.0;
  } else if (r > cut_hi) {
    return 2 * eps * b_hi * (1 - cut_hc / r);
  } else if (r > cut_lo) {
    KK_FLOAT tmp = Kokkos::expf(-(r - cut_0) * a);
    return 2 * eps * (1 - tmp) * tmp * a / r;
  } else if (r > cut_lc) {
    return 2 * eps * b_lo * (1 - cut_lc / r);
  } else {
    return 0.0;
  }
}

/* ----------------------------------------------------------------------
   f2 modulation factor
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT F2_KK(KK_FLOAT r, KK_FLOAT k, KK_FLOAT cut_0, KK_FLOAT cut_lc,
                     KK_FLOAT cut_hc, KK_FLOAT cut_lo, KK_FLOAT cut_hi,
                     KK_FLOAT b_lo, KK_FLOAT b_hi, KK_FLOAT cut_c)
{
  if (r < cut_lc || r > cut_hc) {
    return 0;
  } else if (r < cut_lo) {
    return k * b_lo * (cut_lc - r) * (cut_lc - r);
  } else if (r < cut_hi) {
    return k * 0.5 * ((r - cut_0) * (r - cut_0) - (cut_0 - cut_c) * (cut_0 - cut_c));
  } else {
    return k * b_hi * (cut_hc - r) * (cut_hc - r);
  }
}

KOKKOS_INLINE_FUNCTION
static KK_FLOAT F2_KK(KK_FLOAT r, KK_FLOAT k, KK_FLOAT cut_0, KK_FLOAT cut_lc,
                     KK_FLOAT cut_hc, KK_FLOAT cut_lo, KK_FLOAT cut_hi,
                     KK_FLOAT b_lo, KK_FLOAT b_hi, KK_FLOAT cut_c, KK_FLOAT &df2)
{
  if (r < cut_lc || r > cut_hc) {
    df2 = 0.0;
    return 0.0;
  } else if (r < cut_lo) {
    df2 = 2 * k * b_lo * (r - cut_lc);
    return k * b_lo * (cut_lc - r) * (cut_lc - r);
  } else if (r < cut_hi) {
    df2 = k * (r - cut_0);
    return k * 0.5 * Kokkos::fma((r - cut_0), (r - cut_0), -(cut_0 - cut_c) * (cut_0 - cut_c));
  } else {
    df2 = 2 * k * b_hi * (r - cut_hc);
    return k * b_hi * (cut_hc - r) * (cut_hc - r);
  }
}

/* ----------------------------------------------------------------------
   derivative of f2 modulation factor
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT DF2_KK(KK_FLOAT r, KK_FLOAT k, KK_FLOAT cut_0, KK_FLOAT cut_lc,
                      KK_FLOAT cut_hc, KK_FLOAT cut_lo, KK_FLOAT cut_hi,
                      KK_FLOAT b_lo, KK_FLOAT b_hi)
{
  if (r < cut_lc || r > cut_hc) {
    return 0;
  } else if (r < cut_lo) {
    return 2 * k * b_lo * (r - cut_lc);
  } else if (r < cut_hi) {
    return k * (r - cut_0);
  } else {
    return 2 * k * b_hi * (r - cut_hc);
  }
}

/* ----------------------------------------------------------------------
   f3 modulation factor
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT F3_KK(KK_FLOAT rsq, KK_FLOAT cutsq_ast, KK_FLOAT cut_c,
                     KK_FLOAT lj1, KK_FLOAT lj2, KK_FLOAT eps, KK_FLOAT b,
                     KK_ACC_FLOAT &fpair)
{
  KK_FLOAT evdwl = 0.0;

  if (rsq < cutsq_ast) {
    KK_FLOAT r2inv = 1.0 / rsq;
    KK_FLOAT r6inv = r2inv * r2inv * r2inv;
    fpair = r2inv * r6inv * (12 * lj1 * r6inv - 6 * lj2);
    evdwl = r6inv * (lj1 * r6inv - lj2);
  } else {
    KK_FLOAT r = Kokkos::sqrtf(rsq);
    KK_FLOAT rinv = 1.0 / r;
    fpair = 2 * eps * b * (cut_c * rinv - 1);
    evdwl = eps * b * (cut_c - r) * (cut_c - r);
  }
  return evdwl;
}

/* ----------------------------------------------------------------------
   f4 modulation factor
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT F4_KK(KK_FLOAT theta, KK_FLOAT a, KK_FLOAT theta_0,
                     KK_FLOAT dtheta_ast, KK_FLOAT b, KK_FLOAT dtheta_c)
{
  KK_FLOAT dtheta = theta - theta_0;

  if (Kokkos::fabs(dtheta) > dtheta_c) {
    return 0.0;
  } else if (dtheta > dtheta_ast) {
    return b * (dtheta - dtheta_c) * (dtheta - dtheta_c);
  } else if (dtheta > -dtheta_ast) {
    return 1 - a * dtheta * dtheta;
  } else {
    return b * (dtheta + dtheta_c) * (dtheta + dtheta_c);
  }
}

KOKKOS_INLINE_FUNCTION
static KK_FLOAT F4_KK(KK_FLOAT theta, KK_FLOAT a, KK_FLOAT theta_0,
                     KK_FLOAT dtheta_ast, KK_FLOAT b, KK_FLOAT dtheta_c, KK_FLOAT &df4)
{
  KK_FLOAT dtheta = theta - theta_0;

  if (Kokkos::fabs(dtheta) > dtheta_c) {
    df4 = 0.0;
    return 0.0;
  } else if (dtheta > dtheta_ast) {
    df4 = 2 * b * (dtheta - dtheta_c);
    return b * (dtheta - dtheta_c) * (dtheta - dtheta_c);
  } else if (dtheta > -dtheta_ast) {
    df4 = -2 * a * dtheta;
    return 1 - a * dtheta * dtheta;
  } else {
    df4 = 2 * b * (dtheta + dtheta_c);
    return b * (dtheta + dtheta_c) * (dtheta + dtheta_c);
  }
}

/* ----------------------------------------------------------------------
   derivative of f4 modulation factor

   NOTE: We handle the sin(theta) factor from the partial derivative
   of d(cos(theta))/dtheta externally. The reason for this is
   because the sign of DF4 depends on the sign of theta in the
   function call. It is also more efficient to store sin(theta).
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT DF4_KK(KK_FLOAT theta, KK_FLOAT a, KK_FLOAT theta_0,
                      KK_FLOAT dtheta_ast, KK_FLOAT b, KK_FLOAT dtheta_c)
{
  KK_FLOAT dtheta = theta - theta_0;

  if (Kokkos::fabs(dtheta) > dtheta_c) {
    return 0.0;
  } else if (dtheta > dtheta_ast) {
    return 2 * b * (dtheta - dtheta_c);
  } else if (dtheta > -dtheta_ast) {
    return -2 * a * dtheta;
  } else {
    return 2 * b * (dtheta + dtheta_c);
  }
}

/* ----------------------------------------------------------------------
   f5 modulation factor
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT F5_KK(KK_FLOAT x, KK_FLOAT a, KK_FLOAT x_ast,
                     KK_FLOAT b, KK_FLOAT x_c)
{
  if (x >= 0) {
    return 1.0;
  } else if (x > x_ast) {
    return 1 - a * x * x;
  } else if (x > x_c) {
    return b * (x - x_c) * (x - x_c);
  } else {
    return 0.0;
  }
}

/* ----------------------------------------------------------------------
   derivative of f5 modulation factor
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT DF5_KK(KK_FLOAT x, KK_FLOAT a, KK_FLOAT x_ast,
                      KK_FLOAT b, KK_FLOAT x_c)
{
  if (x >= 0) {
    return 0.0;
  } else if (x > x_ast) {
    return -2 * a * x;
  } else if (x > x_c) {
    return 2 * b * (x - x_c);
  } else {
    return 0.0;
  }
}

/* ----------------------------------------------------------------------
   f6 modulation factor
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT F6_KK(KK_FLOAT theta, KK_FLOAT a, KK_FLOAT b)
{
  if (theta < b) {
    return 0.0;
  } else {
    return 0.5 * a * (theta - b) * (theta - b);
  }
}

KOKKOS_INLINE_FUNCTION
static KK_FLOAT F6_KK(KK_FLOAT theta, KK_FLOAT a, KK_FLOAT b, KK_FLOAT &df6)
{
  if (theta < b) {
    df6 = 0.0;
    return 0.0;
  } else {
    df6 = a * (theta - b);
    return 0.5 * a * (theta - b) * (theta - b);
  }
}

/* ----------------------------------------------------------------------
   derivative of f6 modulation factor
   ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
static KK_FLOAT DF6_KK(KK_FLOAT theta, KK_FLOAT a, KK_FLOAT b)
{
  if (theta < b) {
    return 0.0;
  } else {
    return a * (theta - b);
  }
}

}    // namespace MFOxdnaKokkos

#endif

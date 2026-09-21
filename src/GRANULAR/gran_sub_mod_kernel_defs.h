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

/* ----------------------------------------------------------------------
   Shared definitions for the granular sub-model kernels.

   The gran_sub_mod_*_kernel.h headers hold the arithmetic of the granular
   sub-models as free functions with no state of their own, so that the
   host classes in gran_sub_mod_*.cpp and the accelerator variants of
   pair granular evaluate the exact same expressions.  Only coefficient
   parsing, mixing and restart handling stay in the host classes.

   The headers deliberately depend on nothing but <cmath> and math_const.h
   so that they can be included from device code.  A file that wants the
   device-callable form includes a Kokkos header (e.g. kokkos_type.h)
   *before* this one; GRAN_KERNEL_FN then expands to
   KOKKOS_INLINE_FUNCTION instead of plain inline.

   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#ifndef LMP_GRAN_SUB_MOD_KERNEL_DEFS_H
#define LMP_GRAN_SUB_MOD_KERNEL_DEFS_H

#include "math_const.h"

#ifdef KOKKOS_INLINE_FUNCTION
#define GRAN_KERNEL_FN KOKKOS_INLINE_FUNCTION
#define GRAN_MATH Kokkos
#else
#include <cmath>
#define GRAN_KERNEL_FN inline
#define GRAN_MATH std
#endif

namespace LAMMPS_NS::Granular_NS {

// sub-model identifiers used by the accelerator kernels to select a model
// at run time.  The host classes call the per-model functions directly and
// do not need these.  Keep in sync with gran_sub_mod_register.cpp.

enum GranNormalModel {
  GRAN_NORMAL_NONE = 0,
  GRAN_NORMAL_HOOKE,
  GRAN_NORMAL_HERTZ,
  GRAN_NORMAL_HERTZ_MATERIAL,
  GRAN_NORMAL_DMT,
  GRAN_NORMAL_JKR,
  GRAN_NORMAL_MDR
};

enum GranDampingModel {
  GRAN_DAMPING_NONE = 0,
  GRAN_DAMPING_VELOCITY,
  GRAN_DAMPING_MASS_VELOCITY,
  GRAN_DAMPING_VISCOELASTIC,
  GRAN_DAMPING_TSUJI,
  GRAN_DAMPING_MDR
};

enum GranTangentialModel {
  GRAN_TANGENTIAL_NONE = 0,
  GRAN_TANGENTIAL_LINEAR_NOHISTORY,
  GRAN_TANGENTIAL_LINEAR_HISTORY,
  GRAN_TANGENTIAL_LINEAR_HISTORY_CLASSIC,
  GRAN_TANGENTIAL_MINDLIN_CLASSIC,
  GRAN_TANGENTIAL_MINDLIN
};

enum GranRollingModel {
  GRAN_ROLLING_NONE = 0,
  GRAN_ROLLING_SDS
};

enum GranTwistingModel {
  GRAN_TWISTING_NONE = 0,
  GRAN_TWISTING_MARSHALL,
  GRAN_TWISTING_SDS
};

enum GranHeatModel {
  GRAN_HEAT_NONE = 0,
  GRAN_HEAT_RADIUS,
  GRAN_HEAT_AREA
};

// tolerance below which a history vector is considered to lie in the
// tangential plane already.  Matches EPSILON in gran_sub_mod_tangential.cpp
// and gran_sub_mod_rolling.cpp.

static constexpr double GRAN_KERNEL_EPSILON = 1e-10;

namespace GranKernel {

  // small helpers; local copies so the headers stay free of LAMMPS includes

  template <class T> GRAN_KERNEL_FN T gk_max(const T a, const T b)
  {
    return (a > b) ? a : b;
  }

  template <class T> GRAN_KERNEL_FN T gk_min(const T a, const T b)
  {
    return (a < b) ? a : b;
  }

  template <class T> GRAN_KERNEL_FN T gk_square(const T a)
  {
    return a * a;
  }

  template <class T> GRAN_KERNEL_FN T gk_cube(const T a)
  {
    return a * a * a;
  }

  template <class T> GRAN_KERNEL_FN T gk_dot3(const T *a, const T *b)
  {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  template <class T> GRAN_KERNEL_FN T gk_len3(const T *a)
  {
    return GRAN_MATH::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  }

  template <class T> GRAN_KERNEL_FN void gk_zero3(T *a)
  {
    a[0] = a[1] = a[2] = (T) 0.0;
  }

  template <class T> GRAN_KERNEL_FN void gk_copy3(const T *a, T *b)
  {
    b[0] = a[0];
    b[1] = a[1];
    b[2] = a[2];
  }

  template <class T> GRAN_KERNEL_FN void gk_scale3(const T s, T *a)
  {
    a[0] *= s;
    a[1] *= s;
    a[2] *= s;
  }

  template <class T> GRAN_KERNEL_FN void gk_scale3(const T s, const T *a, T *b)
  {
    b[0] = s * a[0];
    b[1] = s * a[1];
    b[2] = s * a[2];
  }

  template <class T> GRAN_KERNEL_FN void gk_add3(const T *a, const T *b, T *c)
  {
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
  }

  template <class T> GRAN_KERNEL_FN void gk_sub3(const T *a, const T *b, T *c)
  {
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
  }

  template <class T>
  GRAN_KERNEL_FN void gk_scaleadd3(const T s1, const T *a, const T s2, const T *b, T *c)
  {
    c[0] = s1 * a[0] + s2 * b[0];
    c[1] = s1 * a[1] + s2 * b[1];
    c[2] = s1 * a[2] + s2 * b[2];
  }

  template <class T> GRAN_KERNEL_FN void gk_cross3(const T *a, const T *b, T *c)
  {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
  }

  /* --------------------------------------------------------------------
     project v into the plane normal to n, then rescale it back to its
     original magnitude.  Shared by the tangential and rolling models;
     mirrors GranSubMod::rotate_rescale_vec().  Note that a vector that
     lies exactly along n collapses to zero, which is the historical
     behavior and must be preserved.
     -------------------------------------------------------------------- */

  template <class T> GRAN_KERNEL_FN void gk_rotate_rescale_vec(T *v, const T *n)
  {
    const T rsht = gk_dot3(v, n);
    const T shrmag = gk_len3(v);

    v[0] -= rsht * n[0];
    v[1] -= rsht * n[1];
    v[2] -= rsht * n[2];

    const T prjmag = gk_len3(v);
    const T scale = (prjmag > (T) 0.0) ? shrmag / prjmag : (T) 0.0;
    gk_scale3(scale, v);
  }

}    // namespace GranKernel

}    // namespace LAMMPS_NS::Granular_NS

#endif

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
   Stateless kernels for the granular normal sub-models.  See
   gran_sub_mod_kernel_defs.h for how these headers are meant to be used.

   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#ifndef LMP_GRAN_SUB_MOD_NORMAL_KERNEL_H
#define LMP_GRAN_SUB_MOD_NORMAL_KERNEL_H

#include "gran_sub_mod_kernel_defs.h"

namespace LAMMPS_NS::Granular_NS::GranKernel {

static constexpr double GK_PI27SQ = 266.479318829412648029;    // 27*PI^2
static constexpr double GK_THREEROOT3 = 5.19615242270663202362;    // 3*sqrt(3)
static constexpr double GK_SIXROOT6 = 14.69693845669906728801;     // 6*sqrt(6)
static constexpr double GK_INVROOT6 = 0.40824829046386307274;      // 1/sqrt(6)
static constexpr double GK_JKRPREFIX = 1.2277228507842888;         // cbrt(3*PI^2/16)

// coefficients of a normal sub-model that stay constant over a run.
// Emix is only meaningful for the cohesive models.

template <class T> struct GranNormalParams {
  int model;
  T k;
  T cohesion;
  T Emix;
};

/* ----------------------------------------------------------------------
   contact radius
------------------------------------------------------------------------- */

// GranSubModNormal::calculate_contact_radius()

template <class T> GRAN_KERNEL_FN T gran_normal_contact_radius_default(const T dR)
{
  return GRAN_MATH::sqrt(dR);
}

// GranSubModNormalJKR::calculate_contact_radius()

template <class T>
GRAN_KERNEL_FN T gran_normal_jkr_contact_radius(const T Reff, const T dR, const T cohesion,
                                                const T Emix)
{
  const T R2 = Reff * Reff;
  const T dR2 = dR * dR;
  const T t0 = cohesion * cohesion * R2 * R2 * Emix;
  const T t1 = (T) GK_PI27SQ * t0;
  const T t2 = (T) 8.0 * dR * dR2 * Emix * Emix * Emix;
  const T t3 = (T) 4.0 * dR2 * Emix;

  // in case sqrt(0) < 0 due to precision issues
  const T sqrt1 = gk_max((T) 0.0, t0 * (t1 + (T) 2.0 * t2));
  const T t4 = GRAN_MATH::cbrt(t1 + t2 + (T) GK_THREEROOT3 * (T) MathConst::MY_PI *
                               GRAN_MATH::sqrt(sqrt1));
  const T t5 = t3 / t4 + t4 / Emix;
  const T sqrt2 = gk_max((T) 0.0, (T) 2.0 * dR + t5);
  const T t6 = GRAN_MATH::sqrt(sqrt2);
  const T sqrt3 = gk_max((T) 0.0, (T) 4.0 * dR - t5 +
                         (T) GK_SIXROOT6 * cohesion * (T) MathConst::MY_PI * R2 / (Emix * t6));

  return (T) GK_INVROOT6 * (t6 + GRAN_MATH::sqrt(sqrt3));
}

template <class T>
GRAN_KERNEL_FN T gran_normal_contact_radius(const GranNormalParams<T> &p, const T Reff, const T dR)
{
  if (p.model == GRAN_NORMAL_JKR)
    return gran_normal_jkr_contact_radius(Reff, dR, p.cohesion, p.Emix);
  return gran_normal_contact_radius_default(dR);
}

/* ----------------------------------------------------------------------
   touch test
------------------------------------------------------------------------- */

// GranSubModNormalJKR::pulloff_distance()

template <class T>
GRAN_KERNEL_FN T gran_normal_jkr_pulloff_distance(const T Reff, const T cohesion, const T Emix)
{
  if (Reff <= (T) 0.0) return (T) 0.0;
  // defined as positive so center-to-center separation is > radsum
  return (T) GK_JKRPREFIX * GRAN_MATH::cbrt(Reff * cohesion * cohesion / (Emix * Emix));
}

// GranSubModNormalJKR::touch(); touch_prev is the stored contact flag

template <class T>
GRAN_KERNEL_FN bool gran_normal_jkr_touch(const T rsq, const T radsum, const T Reff,
                                          const T cohesion, const T Emix, const bool touch_prev)
{
  if (touch_prev) {
    const T dist_pulloff = radsum + gran_normal_jkr_pulloff_distance(Reff, cohesion, Emix);
    return rsq < (dist_pulloff * dist_pulloff);
  }
  return rsq < (radsum * radsum);
}

template <class T>
GRAN_KERNEL_FN bool gran_normal_touch(const GranNormalParams<T> &p, const T rsq, const T radsum,
                                      const T Reff, const bool touch_prev)
{
  if (p.model == GRAN_NORMAL_JKR)
    return gran_normal_jkr_touch(rsq, radsum, Reff, p.cohesion, p.Emix, touch_prev);
  return rsq < (radsum * radsum);
}

/* ----------------------------------------------------------------------
   normal forces
------------------------------------------------------------------------- */

// GranSubModNormalHooke::calculate_forces()

template <class T> GRAN_KERNEL_FN T gran_normal_hooke_force(const T k, const T delta)
{
  return k * delta;
}

// GranSubModNormalHertz::calculate_forces(); also hertz/material

template <class T>
GRAN_KERNEL_FN T gran_normal_hertz_force(const T k, const T contact_radius, const T delta)
{
  return k * contact_radius * delta;
}

// GranSubModNormalDMT::calculate_forces(); Fne and F_pulloff are kept for set_fncrit()

template <class T>
GRAN_KERNEL_FN T gran_normal_dmt_force(const T k, const T contact_radius, const T delta,
                                       const T cohesion, const T Reff, T &Fne, T &F_pulloff)
{
  Fne = k * contact_radius * delta;
  F_pulloff = (T) 4.0 * (T) MathConst::MY_PI * cohesion * Reff;
  Fne -= F_pulloff;
  return Fne;
}

// GranSubModNormalJKR::calculate_forces()

template <class T>
GRAN_KERNEL_FN T gran_normal_jkr_force(const T k, const T contact_radius, const T Reff,
                                       const T cohesion, const T Emix, T &Fne, T &F_pulloff)
{
  const T a2 = contact_radius * contact_radius;
  Fne = k * contact_radius * a2 / Reff -
      (T) MathConst::MY_2PI * a2 *
          GRAN_MATH::sqrt((T) 4.0 * cohesion * Emix / ((T) MathConst::MY_PI * contact_radius));
  F_pulloff = (T) 3.0 * (T) MathConst::MY_PI * cohesion * Reff;
  return Fne;
}

template <class T>
GRAN_KERNEL_FN T gran_normal_force(const GranNormalParams<T> &p, const T contact_radius,
                                   const T delta, const T Reff, T &Fne, T &F_pulloff)
{
  Fne = (T) 0.0;
  F_pulloff = (T) 0.0;

  switch (p.model) {
    case GRAN_NORMAL_HOOKE:
      return gran_normal_hooke_force(p.k, delta);
    case GRAN_NORMAL_HERTZ:
    case GRAN_NORMAL_HERTZ_MATERIAL:
      return gran_normal_hertz_force(p.k, contact_radius, delta);
    case GRAN_NORMAL_DMT:
      return gran_normal_dmt_force(p.k, contact_radius, delta, p.cohesion, Reff, Fne, F_pulloff);
    case GRAN_NORMAL_JKR:
      return gran_normal_jkr_force(p.k, contact_radius, Reff, p.cohesion, p.Emix, Fne, F_pulloff);
    default:
      return (T) 0.0;
  }
}

/* ----------------------------------------------------------------------
   critical normal force used by the tangential, rolling and twisting models
------------------------------------------------------------------------- */

// GranSubModNormal::set_fncrit()

template <class T> GRAN_KERNEL_FN T gran_normal_fncrit_default(const T Fntot)
{
  return GRAN_MATH::fabs(Fntot);
}

// GranSubModNormalDMT::set_fncrit() and GranSubModNormalJKR::set_fncrit()

template <class T> GRAN_KERNEL_FN T gran_normal_fncrit_cohesive(const T Fne, const T F_pulloff)
{
  return GRAN_MATH::fabs(Fne + (T) 2.0 * F_pulloff);
}

template <class T>
GRAN_KERNEL_FN T gran_normal_fncrit(const GranNormalParams<T> &p, const T Fntot, const T Fne,
                                    const T F_pulloff)
{
  if ((p.model == GRAN_NORMAL_DMT) || (p.model == GRAN_NORMAL_JKR))
    return gran_normal_fncrit_cohesive(Fne, F_pulloff);
  return gran_normal_fncrit_default(Fntot);
}

}    // namespace LAMMPS_NS::Granular_NS::GranKernel

#endif

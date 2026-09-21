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
   Stateless kernels for the granular damping sub-models.  See
   gran_sub_mod_kernel_defs.h for how these headers are meant to be used.

   Every damping model produces a prefactor; the damping force is always
   -prefactor * vnnr.  The prefactor is also read back by the tangential,
   rolling and twisting models, so it is returned rather than folded in.

   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#ifndef LMP_GRAN_SUB_MOD_DAMPING_KERNEL_H
#define LMP_GRAN_SUB_MOD_DAMPING_KERNEL_H

#include "gran_sub_mod_kernel_defs.h"

namespace LAMMPS_NS::Granular_NS::GranKernel {

template <class T> struct GranDampingParams {
  int model;
  T damp;
};

// GranSubModDampingVelocity::calculate_forces()

template <class T> GRAN_KERNEL_FN T gran_damping_velocity_prefactor(const T damp)
{
  return damp;
}

// GranSubModDampingMassVelocity::calculate_forces()

template <class T>
GRAN_KERNEL_FN T gran_damping_mass_velocity_prefactor(const T damp, const T meff)
{
  return damp * meff;
}

// GranSubModDampingViscoelastic::calculate_forces()

template <class T>
GRAN_KERNEL_FN T gran_damping_viscoelastic_prefactor(const T damp, const T meff,
                                                     const T contact_radius)
{
  return damp * meff * contact_radius;
}

// GranSubModDampingTsuji::calculate_forces(); also coeff_restitution, which
// differs only in how damp is derived at setup time

template <class T>
GRAN_KERNEL_FN T gran_damping_tsuji_prefactor(const T damp, const T meff, const T Fnormal,
                                              const T delta)
{
  // in case argument <= 0 due to precision issues
  T sqrt1;
  if (delta > (T) 0.0)
    sqrt1 = gk_max((T) 0.0, meff * Fnormal / delta);
  else
    sqrt1 = (T) 0.0;
  return damp * GRAN_MATH::sqrt(sqrt1);
}

template <class T>
GRAN_KERNEL_FN T gran_damping_prefactor(const GranDampingParams<T> &p, const T meff,
                                        const T contact_radius, const T Fnormal, const T delta)
{
  switch (p.model) {
    case GRAN_DAMPING_VELOCITY:
      return gran_damping_velocity_prefactor(p.damp);
    case GRAN_DAMPING_MASS_VELOCITY:
      return gran_damping_mass_velocity_prefactor(p.damp, meff);
    case GRAN_DAMPING_VISCOELASTIC:
      return gran_damping_viscoelastic_prefactor(p.damp, meff, contact_radius);
    case GRAN_DAMPING_TSUJI:
      return gran_damping_tsuji_prefactor(p.damp, meff, Fnormal, delta);
    default:
      return (T) 0.0;
  }
}

}    // namespace LAMMPS_NS::Granular_NS::GranKernel

#endif

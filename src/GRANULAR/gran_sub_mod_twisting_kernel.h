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
   Stateless kernels for the granular twisting sub-models.  See
   gran_sub_mod_kernel_defs.h for how these headers are meant to be used.

   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#ifndef LMP_GRAN_SUB_MOD_TWISTING_KERNEL_H
#define LMP_GRAN_SUB_MOD_TWISTING_KERNEL_H

#include "gran_sub_mod_kernel_defs.h"

namespace LAMMPS_NS::Granular_NS::GranKernel {

template <class T> struct GranTwistingParams {
  int model;
  T k;          // sds only
  T damp;       // sds only
  T mu;         // sds only
  T k_tang;     // marshall only, from the tangential model
  T mu_tang;    // marshall only, from the tangential model
};

/* ----------------------------------------------------------------------
   shared core of both twisting models, eqs. 30, 34 and 44 of Marshall
------------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN T gran_twisting_core(const T k, const T damp, const T mu, const T magtwist,
                                    const T dt, const T Fncrit, const int history_update,
                                    T *history)
{
  if (history_update) history[0] += magtwist * dt;

  // M_t torque (eq 30)
  T magtortwist = -k * history[0] - damp * magtwist;
  const T signtwist = (T) ((magtwist > (T) 0.0) - (magtwist < (T) 0.0));
  const T Mtcrit = mu * Fncrit;    // critical torque (eq 44)

  if (GRAN_MATH::fabs(magtortwist) > Mtcrit) {
    history[0] = (Mtcrit * signtwist - damp * magtwist) / k;
    magtortwist = -Mtcrit * signtwist;    // eq 34
  }

  return magtortwist;
}

/* ----------------------------------------------------------------------
   GranSubModTwistingMarshall::calculate_forces(); the twist coefficients
   follow from the tangential model and the contact geometry (eq 32)
------------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN T gran_twisting_marshall(const GranTwistingParams<T> &p, const T magtwist,
                                        const T dt, const T contact_radius, const T Fncrit,
                                        const T tangential_damp, const int history_update,
                                        T *history)
{
  const T k = (T) 0.5 * p.k_tang * contact_radius * contact_radius;
  const T damp = (T) 0.5 * tangential_damp * contact_radius * contact_radius;
  const T mu = (T) MathConst::TWOTHIRDS * p.mu_tang * contact_radius;

  return gran_twisting_core(k, damp, mu, magtwist, dt, Fncrit, history_update, history);
}

/* ----------------------------------------------------------------------
   GranSubModTwistingSDS::calculate_forces()
------------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN T gran_twisting_sds(const GranTwistingParams<T> &p, const T magtwist, const T dt,
                                   const T Fncrit, const int history_update, T *history)
{
  return gran_twisting_core(p.k, p.damp, p.mu, magtwist, dt, Fncrit, history_update, history);
}

/* ---------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN T gran_twisting_forces(const GranTwistingParams<T> &p, const T magtwist, const T dt,
                                      const T contact_radius, const T Fncrit,
                                      const T tangential_damp, const int history_update,
                                      T *history)
{
  switch (p.model) {
    case GRAN_TWISTING_MARSHALL:
      return gran_twisting_marshall(p, magtwist, dt, contact_radius, Fncrit, tangential_damp,
                                    history_update, history);
    case GRAN_TWISTING_SDS:
      return gran_twisting_sds(p, magtwist, dt, Fncrit, history_update, history);
    default:
      return (T) 0.0;
  }
}

}    // namespace LAMMPS_NS::Granular_NS::GranKernel

#endif

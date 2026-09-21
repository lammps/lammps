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
   Stateless kernels for the granular rolling sub-models.  See
   gran_sub_mod_kernel_defs.h for how these headers are meant to be used.

   Unlike the tangential models, the rolling model works on a copy of the
   history and only stores it back when history_update is set.  That
   difference is preserved here.

   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#ifndef LMP_GRAN_SUB_MOD_ROLLING_KERNEL_H
#define LMP_GRAN_SUB_MOD_ROLLING_KERNEL_H

#include "gran_sub_mod_kernel_defs.h"

namespace LAMMPS_NS::Granular_NS::GranKernel {

template <class T> struct GranRollingParams {
  int model;
  T k;
  T gamma;
  T mu;
};

template <class T> struct GranRollingState {
  const T *nx;
  const T *nx_unrotated;
  const T *vrl;
  T dt;
  T Fncrit;
  int synchronized_verlet;
  int history_update;
};

/* ----------------------------------------------------------------------
   GranSubModRollingSDS::calculate_forces()
------------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN void gran_rolling_sds(const GranRollingParams<T> &p, const GranRollingState<T> &s,
                                     T *history, T *fr)
{
  int frameupdate;
  T hist_temp[3], temp_array[3];

  const T *nx = s.nx;
  const T *vrl = s.vrl;
  const T k = p.k;

  const T Frcrit = p.mu * s.Fncrit;

  gk_copy3(history, hist_temp);

  if (s.history_update) {
    T rolldotn = gk_dot3(hist_temp, nx);

    frameupdate = (GRAN_MATH::fabs(rolldotn) * k) > ((T) GRAN_KERNEL_EPSILON * Frcrit);
    if (frameupdate) gk_rotate_rescale_vec(hist_temp, nx);

    // update history at half-step
    gk_scale3(s.dt, vrl, temp_array);
    gk_add3(hist_temp, temp_array, hist_temp);

    // rotate into tangential plane at full-step for synchronized_verlet
    if (s.synchronized_verlet == 1) {
      rolldotn = gk_dot3(hist_temp, s.nx_unrotated);
      frameupdate = (GRAN_MATH::fabs(rolldotn) * k) > ((T) GRAN_KERNEL_EPSILON * Frcrit);
      if (frameupdate) gk_rotate_rescale_vec(hist_temp, s.nx_unrotated);
    }
  }

  gk_scaleadd3(-k, hist_temp, -p.gamma, vrl, fr);

  // rescale frictional displacements and forces if needed

  const T magfr = gk_len3(fr);
  if (magfr > Frcrit) {
    const T rollmag = gk_len3(hist_temp);
    if (rollmag != (T) 0.0) {
      const T k_inv = (T) 1.0 / k;
      const T magfr_inv = (T) 1.0 / magfr;
      gk_scale3(-Frcrit * k_inv * magfr_inv, fr, hist_temp);
      gk_scale3(-p.gamma * k_inv, vrl, temp_array);
      gk_add3(hist_temp, temp_array, hist_temp);

      gk_scale3(Frcrit * magfr_inv, fr);
    } else {
      gk_zero3(fr);
    }
  }

  if (s.history_update) gk_copy3(hist_temp, history);
}

/* ---------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN void gran_rolling_forces(const GranRollingParams<T> &p,
                                        const GranRollingState<T> &s, T *history, T *fr)
{
  if (p.model == GRAN_ROLLING_SDS)
    gran_rolling_sds(p, s, history, fr);
  else
    gk_zero3(fr);
}

}    // namespace LAMMPS_NS::Granular_NS::GranKernel

#endif

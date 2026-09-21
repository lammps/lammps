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
   Stateless kernels for the granular tangential sub-models.  See
   gran_sub_mod_kernel_defs.h for how these headers are meant to be used.

   The history pointer is written in place.  Note that the Coulomb rescale
   at the end of each history model updates the stored history even when
   history_update is 0; this is long-standing behavior and callers must
   not gate the write-back on history_update.

   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#ifndef LMP_GRAN_SUB_MOD_TANGENTIAL_KERNEL_H
#define LMP_GRAN_SUB_MOD_TANGENTIAL_KERNEL_H

#include "gran_sub_mod_kernel_defs.h"

namespace LAMMPS_NS::Granular_NS::GranKernel {

// coefficients of a tangential sub-model that stay constant over a run

template <class T> struct GranTangentialParams {
  int model;
  T k;
  T xt;
  T mu;
  int mindlin_force;
  int mindlin_rescale;
  int contact_radius_flag;
};

// per-contact inputs shared by all tangential sub-models

template <class T> struct GranTangentialState {
  const T *nx;
  const T *nx_unrotated;
  const T *vtr;
  T vrel;
  T dt;
  T contact_radius;
  T damp_prefactor;
  T Fncrit;
  int synchronized_verlet;
  int history_update;
};

// tangential damping coefficient; also read back by the Marshall twisting model

template <class T>
GRAN_KERNEL_FN T gran_tangential_damp(const GranTangentialParams<T> &p, const T damp_prefactor)
{
  return p.xt * damp_prefactor;
}

/* ----------------------------------------------------------------------
   GranSubModTangentialLinearNoHistory::calculate_forces()
------------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN void gran_tangential_linear_nohistory(const GranTangentialParams<T> &p,
                                                     const GranTangentialState<T> &s, T *fs)
{
  const T damp = gran_tangential_damp(p, s.damp_prefactor);
  const T Fscrit = p.mu * s.Fncrit;
  const T fsmag = damp * s.vrel;

  T Ft;
  if (s.vrel != (T) 0.0)
    Ft = gk_min(Fscrit, fsmag) / s.vrel;
  else
    Ft = (T) 0.0;

  gk_scale3(-Ft, s.vtr, fs);
}

/* ----------------------------------------------------------------------
   GranSubModTangentialLinearHistory::calculate_forces()
------------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN void gran_tangential_linear_history(const GranTangentialParams<T> &p,
                                                   const GranTangentialState<T> &s, T *history,
                                                   T *fs)
{
  T hist_increment[3], fdamp[3], vtr2[3];
  int frame_update = 0;

  const T *nx = s.nx;
  const T *vtr = s.vtr;
  const T k = p.k;

  const T damp = gran_tangential_damp(p, s.damp_prefactor);
  const T Fscrit = s.Fncrit * p.mu;

  // rotate and update displacements
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235

  if (s.history_update) {
    T rsht = gk_dot3(history, nx);
    frame_update = (GRAN_MATH::fabs(rsht) * k) > ((T) GRAN_KERNEL_EPSILON * Fscrit);

    if (frame_update) gk_rotate_rescale_vec(history, nx);

    gk_scale3(s.dt, vtr, hist_increment);
    gk_add3(history, hist_increment, history);

    if (s.synchronized_verlet == 1) {
      rsht = gk_dot3(history, s.nx_unrotated);
      frame_update = (GRAN_MATH::fabs(rsht) * k) > ((T) GRAN_KERNEL_EPSILON * Fscrit);
      if (frame_update) gk_rotate_rescale_vec(history, s.nx_unrotated);
    }
  }

  // tangential forces = history + tangential velocity damping

  gk_scale3(-k, history, fs);

  if (frame_update && (s.synchronized_verlet == 1)) {
    gk_copy3(vtr, vtr2);
    gk_rotate_rescale_vec(vtr2, s.nx_unrotated);
  } else {
    gk_copy3(vtr, vtr2);
  }
  gk_scale3(-damp, vtr2, fdamp);
  gk_add3(fs, fdamp, fs);

  // rescale frictional displacements and forces if needed

  const T magfs = gk_len3(fs);
  if (magfs > Fscrit) {
    const T shrmag = gk_len3(history);
    if (shrmag != (T) 0.0) {
      gk_scale3(Fscrit / magfs, fs);
      gk_sub3(fs, fdamp, history);
      gk_scale3((T) -1.0 / k, history);
    } else {
      gk_zero3(fs);
    }
  }
}

/* ----------------------------------------------------------------------
   GranSubModTangentialLinearHistoryClassic::calculate_forces();
   also covers mindlin_classic, which only sets contact_radius_flag
------------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN void gran_tangential_classic(const GranTangentialParams<T> &p,
                                            const GranTangentialState<T> &s, T *history, T *fs)
{
  T hist_increment[3], fdamp[3];

  const T *nx = s.nx;
  const T *vtr = s.vtr;

  const T damp = gran_tangential_damp(p, s.damp_prefactor);
  const T Fscrit = s.Fncrit * p.mu;

  // update history

  if (s.history_update) {
    gk_scale3(s.dt, vtr, hist_increment);
    gk_add3(history, hist_increment, history);
  }

  const T shrmag = gk_len3(history);

  // rotate shear displacements

  if (s.history_update) {
    const T rsht = gk_dot3(history, nx);
    gk_scale3(rsht, nx, hist_increment);
    gk_sub3(history, hist_increment, history);
  }

  // classic model can only set contact_radius_flag through hertz

  const T k_scaled = p.contact_radius_flag ? (p.k * s.contact_radius) : p.k;
  gk_scale3(-k_scaled, history, fs);

  // damping force, note that damp automatically has a factor
  //   of contact radius with hertz (sets viscoelastic damping)
  //   but not with hooke (sets mass_velocity damping)

  gk_scale3(-damp, vtr, fdamp);
  gk_add3(fs, fdamp, fs);

  // rescale frictional displacements and forces if needed

  const T magfs = gk_len3(fs);
  if (magfs > Fscrit) {
    if (shrmag != (T) 0.0) {
      gk_scale3(Fscrit / magfs, fs);
      gk_sub3(fs, fdamp, history);
      gk_scale3((T) -1.0 / k_scaled, history);
    } else {
      gk_zero3(fs);
    }
  }
}

/* ----------------------------------------------------------------------
   GranSubModTangentialMindlin::calculate_forces(); also covers
   mindlin/force, mindlin_rescale and mindlin_rescale/force.
   history[3] is only touched by the rescale variants.
------------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN void gran_tangential_mindlin(const GranTangentialParams<T> &p,
                                            const GranTangentialState<T> &s, T *history, T *fs)
{
  T hist_increment[3], fdamp[3], vtr2[3];
  int frame_update = 0;

  const T *nx = s.nx;
  const T *vtr = s.vtr;
  const T contact_radius = s.contact_radius;

  const T damp = gran_tangential_damp(p, s.damp_prefactor);
  const T Fscrit = s.Fncrit * p.mu;

  const T k_scaled = p.k * contact_radius;

  // on unloading, rescale the shear displacements/force

  if (p.mindlin_rescale)
    if (contact_radius < history[3]) gk_scale3(contact_radius / history[3], history);

  // rotate and update displacements / force
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235

  if (s.history_update) {
    T rsht = gk_dot3(history, nx);
    if (p.mindlin_force)
      frame_update = GRAN_MATH::fabs(rsht) > ((T) GRAN_KERNEL_EPSILON * Fscrit);
    else
      frame_update = (GRAN_MATH::fabs(rsht) * k_scaled) > ((T) GRAN_KERNEL_EPSILON * Fscrit);

    if (frame_update) gk_rotate_rescale_vec(history, nx);

    // update history
    if (p.mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      gk_scale3(-k_scaled * s.dt, vtr, hist_increment);
    } else {
      gk_scale3(s.dt, vtr, hist_increment);
    }
    gk_add3(history, hist_increment, history);

    if (p.mindlin_rescale) history[3] = contact_radius;

    if (s.synchronized_verlet == 1) {
      // second projection to full step normal
      rsht = gk_dot3(history, s.nx_unrotated);
      if (p.mindlin_force)
        frame_update = GRAN_MATH::fabs(rsht) > ((T) GRAN_KERNEL_EPSILON * Fscrit);
      else
        frame_update = (GRAN_MATH::fabs(rsht) * k_scaled) > ((T) GRAN_KERNEL_EPSILON * Fscrit);
      if (frame_update) gk_rotate_rescale_vec(history, s.nx_unrotated);
    }
  }

  // tangential forces = history + tangential velocity damping

  if (!p.mindlin_force) {
    gk_scale3(-k_scaled, history, fs);
  } else {
    gk_copy3(history, fs);
  }

  // rotating vtr for damping term in nx direction

  if (frame_update && s.synchronized_verlet) {
    gk_copy3(vtr, vtr2);
    gk_rotate_rescale_vec(vtr2, s.nx_unrotated);
  } else {
    gk_copy3(vtr, vtr2);
  }
  gk_scale3(-damp, vtr2, fdamp);
  gk_add3(fs, fdamp, fs);

  // rescale frictional displacements and forces if needed

  const T magfs = gk_len3(fs);
  if (magfs > Fscrit) {
    const T shrmag = gk_len3(history);
    if (shrmag != (T) 0.0) {
      // rescale shear force
      gk_scale3(Fscrit / magfs, fs);

      // set shear to elastic component of rescaled force
      //  may have extra factor of k_scaled that is then removed
      gk_sub3(fs, fdamp, history);
      if (!p.mindlin_force) gk_scale3((T) -1.0 / k_scaled, history);
    } else {
      gk_zero3(fs);
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class T>
GRAN_KERNEL_FN void gran_tangential_forces(const GranTangentialParams<T> &p,
                                           const GranTangentialState<T> &s, T *history, T *fs)
{
  switch (p.model) {
    case GRAN_TANGENTIAL_LINEAR_NOHISTORY:
      gran_tangential_linear_nohistory(p, s, fs);
      break;
    case GRAN_TANGENTIAL_LINEAR_HISTORY:
      gran_tangential_linear_history(p, s, history, fs);
      break;
    case GRAN_TANGENTIAL_LINEAR_HISTORY_CLASSIC:
    case GRAN_TANGENTIAL_MINDLIN_CLASSIC:
      gran_tangential_classic(p, s, history, fs);
      break;
    case GRAN_TANGENTIAL_MINDLIN:
      gran_tangential_mindlin(p, s, history, fs);
      break;
    default:
      gk_zero3(fs);
      break;
  }
}

}    // namespace LAMMPS_NS::Granular_NS::GranKernel

#endif

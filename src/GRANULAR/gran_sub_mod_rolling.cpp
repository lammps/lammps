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

#include "gran_sub_mod_rolling.h"

#include "error.h"
#include "gran_sub_mod_normal.h"
#include "granular_model.h"
#include "math_extra.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathExtra;

static constexpr double EPSILON = 1e-10;

/* ----------------------------------------------------------------------
   Default rolling friction model
------------------------------------------------------------------------- */

GranSubModRolling::GranSubModRolling(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp)
{
  allow_synchronization = 0;
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModRollingNone::GranSubModRollingNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModRolling(gm, lmp)
{
  allow_synchronization = 1;
}

/* ----------------------------------------------------------------------
   SDS rolling friction model
------------------------------------------------------------------------- */

GranSubModRollingSDS::GranSubModRollingSDS(GranularModel *gm, LAMMPS *lmp) :
    GranSubModRolling(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
  allow_synchronization = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModRollingSDS::coeffs_to_local()
{
  k = coeffs[0];
  gamma = coeffs[1];
  mu = coeffs[2];

  if (k < 0.0 || mu < 0.0 || gamma < 0.0) error->all(FLERR, "Illegal SDS rolling model");
}

/* ---------------------------------------------------------------------- */

void GranSubModRollingSDS::calculate_forces()
{
  int rhist0, rhist1, rhist2, frameupdate;
  double Frcrit, rolldotn, rollmag, magfr, hist_temp[3], temp_array[3];
  double k_inv, magfr_inv;

  double *nx = gm->nx;
  double *nx_unrotated = gm->nx_unrotated;
  double *vrl = gm->vrl;
  double *fr = gm->fr;
  double dt = gm->dt;
  double *history = gm->history;
  int history_update = gm->history_update;

  double Fncrit = gm->normal_model->get_fncrit();
  Frcrit = mu * Fncrit;

  rhist0 = history_index;
  rhist1 = rhist0 + 1;
  rhist2 = rhist1 + 1;

  Frcrit = mu * gm->normal_model->get_fncrit();

  hist_temp[0] = history[rhist0];
  hist_temp[1] = history[rhist1];
  hist_temp[2] = history[rhist2];

  if (history_update) {
    rolldotn = dot3(hist_temp, nx);

    frameupdate = (fabs(rolldotn) * k) > (EPSILON * Frcrit);
    if (frameupdate) rotate_rescale_vec(hist_temp, nx);

    // update history at half-step
    scale3(dt, vrl, temp_array);
    add3(hist_temp, temp_array, hist_temp);

    // rotate into tangential plane at full-step for synchronized_verlet
    if (gm->synchronized_verlet == 1) {
      rolldotn = dot3(hist_temp, nx_unrotated);
      frameupdate = (fabs(rolldotn) * k) > (EPSILON * Frcrit);
      if (frameupdate) rotate_rescale_vec(hist_temp, nx_unrotated);
    }
  }

  scaleadd3(-k, hist_temp, -gamma, vrl, fr);

  // rescale frictional displacements and forces if needed
  magfr = len3(fr);
  if (magfr > Frcrit) {
    rollmag = len3(hist_temp);
    if (rollmag != 0.0) {
      k_inv = 1.0 / k;
      magfr_inv = 1.0 / magfr;
      scale3(-Frcrit * k_inv * magfr_inv, fr, hist_temp);
      scale3(-gamma * k_inv, vrl, temp_array);
      add3(hist_temp, temp_array, hist_temp);

      scale3(Frcrit * magfr_inv, fr);
    } else {
      zero3(fr);
    }
  }

  if (history_update) {
    history[rhist0] = hist_temp[0];
    history[rhist1] = hist_temp[1];
    history[rhist2] = hist_temp[2];
  }
}

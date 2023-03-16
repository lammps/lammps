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

#include "gran_sub_mod_damping.h"
#include "gran_sub_mod_normal.h"
#include "gran_sub_mod_tangential.h"
#include "granular_model.h"
#include "error.h"
#include "math_const.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathConst;
using namespace MathExtra;

#define EPSILON 1e-10

/* ----------------------------------------------------------------------
   Default model
------------------------------------------------------------------------- */

GranSubModTangential::GranSubModTangential(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp) {}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModTangentialNone::GranSubModTangentialNone(GranularModel *gm, LAMMPS *lmp) : GranSubModTangential(gm, lmp) {}

/* ----------------------------------------------------------------------
   Linear model with no history
------------------------------------------------------------------------- */

GranSubModTangentialLinearNoHistory::GranSubModTangentialLinearNoHistory(GranularModel *gm, LAMMPS *lmp) : GranSubModTangential(gm, lmp)
{
  num_coeffs = 2;
  size_history = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearNoHistory::coeffs_to_local()
{
  k = 0.0; // No tangential stiffness with no history
  xt = coeffs[0];
  mu = coeffs[1];

  if (k < 0.0 || xt < 0.0 || mu < 0.0)
    error->all(FLERR, "Illegal linear no history tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearNoHistory::calculate_forces()
{
  // classic pair gran/hooke (no history)
  damp = xt * gm->damping_model->damp_prefactor;

  double Fscrit = mu * gm->normal_model->Fncrit;
  double fsmag = damp * gm->vrel;

  double Ft;
  if (gm->vrel != 0.0) Ft = MIN(Fscrit, fsmag) / gm->vrel;
  else Ft = 0.0;

  scale3(-Ft, gm->vtr, gm->fs);
}

/* ----------------------------------------------------------------------
   Linear model with history
------------------------------------------------------------------------- */

GranSubModTangentialLinearHistory::GranSubModTangentialLinearHistory(GranularModel *gm, LAMMPS *lmp) : GranSubModTangential(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistory::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu = coeffs[2];

  if (k < 0.0 || xt < 0.0 || mu < 0.0)
    error->all(FLERR, "Illegal linear tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistory::calculate_forces()
{
  // Note: this is the same as the base Mindlin calculation except k isn't scaled by area
  double magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl, temp_array[3];
  int frame_update = 0;

  damp = xt * gm->damping_model->damp_prefactor;

  double Fscrit = gm->normal_model->Fncrit * mu;
  double *history = & gm->history[history_index];

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    frame_update = (fabs(rsht) * k) > (EPSILON * Fscrit);

    if (frame_update) {
      shrmag = len3(history);

      // projection
      scale3(rsht, gm->nx, temp_array);
      sub3(history, temp_array, history);

      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history, tangential force
    // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
    scale3(gm->dt, gm->vtr, temp_array);
    add3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  scale3(-k, history, gm->fs);
  scale3(damp, gm->vtr, temp_array);
  sub3(gm->fs, temp_array, gm->fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, gm->fs, history);
      scale3(damp, gm->vtr, temp_array);
      add3(history, temp_array, history);
      scale3(-1.0 / k, history);
      scale3(Fscrit * magfs_inv, gm->fs);
    } else {
      zero3(gm->fs);
    }
  }
}

/* ----------------------------------------------------------------------
   Linear model with history from pair gran/hooke/history
------------------------------------------------------------------------- */

GranSubModTangentialLinearHistoryClassic::GranSubModTangentialLinearHistoryClassic(GranularModel *gm, LAMMPS *lmp) : GranSubModTangentialLinearHistory(gm, lmp)
{
  area_flag = 1; // Sets gran/hooke/history behavior
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistoryClassic::calculate_forces()
{
  double magfs, magfs_inv, rsht, shrmag;
  double temp_array[3];

  damp = xt * gm->damping_model->damp_prefactor;

  double Fscrit = gm->normal_model->Fncrit * mu;
  double *history = & gm->history[history_index];

  // update history
  if (gm->history_update) {
    scale3(gm->dt, gm->vtr, temp_array);
    add3(history, temp_array, history);
  }

  shrmag = len3(history);

  // rotate shear displacements
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    scale3(rsht, gm->nx, temp_array);
    sub3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  if (area_flag) scale3(-k * gm->area, history, gm->fs);
  else scale3(-k, history, gm->fs);
  scale3(damp, gm->vtr, temp_array);
  sub3(gm->fs, temp_array, gm->fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  if (magfs > Fscrit) {
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, gm->fs, history);
      scale3(damp, gm->vtr, temp_array);
      add3(history, temp_array, history);
      scale3(-1.0 / k, history);
      scale3(Fscrit * magfs_inv, gm->fs);
    } else {
      zero3(gm->fs);
    }
  }
}

/* ----------------------------------------------------------------------
   Mindlin from pair gran/hertz/history
------------------------------------------------------------------------- */

GranSubModTangentialMindlinClassic::GranSubModTangentialMindlinClassic(GranularModel *gm, LAMMPS *lmp) : GranSubModTangentialLinearHistoryClassic(gm, lmp)
{
  area_flag = 1; // Sets gran/hertz/history behavior
}

/* ----------------------------------------------------------------------
   Mindlin model
------------------------------------------------------------------------- */

GranSubModTangentialMindlin::GranSubModTangentialMindlin(GranularModel *gm, LAMMPS *lmp) : GranSubModTangential(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
  mindlin_force = 0;
  mindlin_rescale = 0;
  area_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu = coeffs[2];

  if (k == -1) {
    if (!gm->normal_model->material_properties)
      error->all(FLERR, "Must either specify tangential stiffness or material properties for normal model for the Mindlin tangential style");

    double Emod = gm->normal_model->Emod;
    double poiss = gm->normal_model->poiss;

    if (gm->contact_type == PAIR) {
      k = 8.0 * mix_stiffnessG(Emod, Emod, poiss, poiss);
    } else {
      k = 8.0 * mix_stiffnessG_wall(Emod, poiss);
    }
  }

  if (k < 0.0 || xt < 0.0 || mu < 0.0)
    error->all(FLERR, "Illegal Mindlin tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::mix_coeffs(double* icoeffs, double* jcoeffs)
{
  if (icoeffs[0] == -1 || jcoeffs[0] == -1) coeffs[0] = -1;
  else coeffs[0] = mix_geom(icoeffs[0], jcoeffs[0]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::calculate_forces()
{
  double k_scaled, magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  damp = xt * gm->damping_model->damp_prefactor;

  double *history = & gm->history[history_index];
  double Fscrit = gm->normal_model->Fncrit * mu;

  k_scaled = k * gm->area;

  // on unloading, rescale the shear displacements/force
  if (mindlin_rescale)
    if (gm->area < history[3])
      scale3(gm->area / history[3], history);

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    if (mindlin_force) {
      frame_update = fabs(rsht) > (EPSILON * Fscrit);
    } else {
      frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
    }

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, gm->nx, temp_array);
      sub3(history, temp_array, history);
      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      scale3(-k_scaled * gm->dt, gm->vtr, temp_array);
    } else {
      scale3(gm->dt, gm->vtr, temp_array);
    }
    add3(history, temp_array, history);

    if (mindlin_rescale) history[3] = gm->area;
  }

  // tangential forces = history + tangential velocity damping
  scale3(-damp, gm->vtr, gm->fs);

  if (!mindlin_force) {
    scale3(k_scaled, history, temp_array);
    sub3(gm->fs, temp_array, gm->fs);
  } else {
    add3(gm->fs, history, gm->fs);
  }

  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, gm->fs, history);
      scale3(damp, gm->vtr, temp_array);
      add3(history, temp_array, history);

      if (!mindlin_force)
        scale3(-1.0 / k_scaled, history);

      scale3(Fscrit * magfs_inv, gm->fs);
    } else {
      zero3(gm->fs);
    }
  }
}

/* ----------------------------------------------------------------------
   Mindlin force model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinForce::GranSubModTangentialMindlinForce(GranularModel *gm, LAMMPS *lmp) : GranSubModTangentialMindlin(gm, lmp)
{
  mindlin_force = 1;
}

/* ----------------------------------------------------------------------
   Mindlin rescale model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinRescale::GranSubModTangentialMindlinRescale(GranularModel *gm, LAMMPS *lmp) : GranSubModTangentialMindlin(gm, lmp)
{
  size_history = 4;
  mindlin_rescale = 1;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

/* ----------------------------------------------------------------------
   Mindlin rescale force model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinRescaleForce::GranSubModTangentialMindlinRescaleForce(GranularModel *gm, LAMMPS *lmp) : GranSubModTangentialMindlin(gm, lmp)
{
  size_history = 4;
  mindlin_force = 1;
  mindlin_rescale = 1;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

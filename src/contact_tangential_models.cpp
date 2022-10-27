/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "contact_damping_models.h"
#include "contact_normal_models.h"
#include "contact_tangential_models.h"
#include "contact.h"
#include "error.h"
#include "math_const.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace Contact;
using namespace MathConst;
using namespace MathExtra;

/* ----------------------------------------------------------------------
   Default model
------------------------------------------------------------------------- */

TangentialModel::TangentialModel(LAMMPS *lmp) : SubModel(lmp) {}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

TangentialNone::TangentialNone(LAMMPS *lmp) : TangentialModel(lmp) {}

/* ----------------------------------------------------------------------
   Linear model with no history
------------------------------------------------------------------------- */

TangentialLinearNoHistory::TangentialLinearNoHistory(LAMMPS *lmp) : TangentialModel(lmp)
{
  num_coeffs = 2;
  size_history = 3;
}

/* ---------------------------------------------------------------------- */

void TangentialLinearNoHistory::coeffs_to_local()
{
  k = 0.0; // No tangential stiffness with no history
  xt = coeffs[0];
  mu = coeffs[1];

  if (k < 0.0 || xt < 0.0 || mu < 0.0)
    error->all(FLERR, "Illegal linear no history tangential model");
}

/* ---------------------------------------------------------------------- */

void TangentialLinearNoHistory::calculate_forces()
{
  // classic pair gran/hooke (no history)
  damp = xt * contact->damping_model->damp_prefactor;

  double Fscrit = mu * contact->normal_model->Fncrit;
  double fsmag = damp * contact->vrel;

  double Ft;
  if (contact->vrel != 0.0) Ft = MIN(Fscrit, fsmag) / contact->vrel;
  else Ft = 0.0;

  scale3(-Ft, contact->vtr, contact->fs);
}

/* ----------------------------------------------------------------------
   Linear model with history
------------------------------------------------------------------------- */

TangentialLinearHistory::TangentialLinearHistory(LAMMPS *lmp) : TangentialModel(lmp)
{
  num_coeffs = 3;
  size_history = 3;
}

/* ---------------------------------------------------------------------- */

void TangentialLinearHistory::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu = coeffs[2];

  if (k < 0.0 || xt < 0.0 || mu < 0.0)
    error->all(FLERR, "Illegal linear tangential model");
}

/* ---------------------------------------------------------------------- */

void TangentialLinearHistory::calculate_forces()
{
  // Note: this is the same as the base Mindlin calculation except k isn't scaled by area
  double magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl, temp_array[3];
  int frame_update = 0;

  damp = xt * contact->damping_model->damp_prefactor;

  double Fscrit = contact->normal_model->Fncrit * mu;
  double *history = & contact->history[history_index];

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (contact->history_update) {
    rsht = dot3(history, contact->nx);
    frame_update = (fabs(rsht) * k) > (EPSILON * Fscrit);

    if (frame_update) {
      shrmag = len3(history);

      // projection
      scale3(rsht, contact->nx, temp_array);
      sub3(history, temp_array, history);

      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history, tangential force
    // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
    scale3(contact->dt, contact->vtr, temp_array);
    add3(history, temp_array, history);
  }


  // tangential forces = history + tangential velocity damping
  scale3(-k, history, contact->fs);
  scale3(damp, contact->vtr, temp_array);
  sub3(contact->fs, temp_array, contact->fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(contact->fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, contact->fs, history);
      scale3(damp, contact->vtr, temp_array);
      add3(history, temp_array, history);
      scale3(-1.0 / k, history);
      scale3(Fscrit * magfs_inv, contact->fs);
    } else {
      zero3(contact->fs);
    }
  }
}

/* ----------------------------------------------------------------------
   Linear model with history from pair gran/hooke/history
------------------------------------------------------------------------- */

TangentialLinearHistoryClassic::TangentialLinearHistoryClassic(LAMMPS *lmp) : TangentialLinearHistory(lmp)
{
  scale_area = 0; // Sets gran/hooke/history behavior
}

/* ---------------------------------------------------------------------- */

void TangentialLinearHistoryClassic::calculate_forces()
{
  double k_scaled, magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  k_scaled = k;
  if (scale_area) k_scaled *= contact->area;

  damp = xt * contact->damping_model->damp_prefactor;

  double Fscrit = contact->normal_model->Fncrit * mu;
  double *history = & contact->history[history_index];

  // update history
  if (contact->history_update) {
    scale3(contact->dt, contact->vtr, temp_array);
    add3(history, temp_array, history);
  }

  shrmag = len3(history);

  // rotate shear displacements
  if (contact->history_update) {
    rsht = dot3(history, contact->nx);
    scale3(rsht, contact->nx, temp_array);
    sub3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  scale3(-k_scaled, history, contact->fs);
  scale3(damp, contact->vtr, temp_array);
  sub3(contact->fs, temp_array, contact->fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(contact->fs);
  if (magfs > Fscrit) {
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, contact->fs, history);
      scale3(damp, contact->vtr, temp_array);
      add3(history, temp_array, history);
      temp_dbl = -1.0 / k_scaled;
      if (scale_area) temp_dbl /= contact->area;
      scale3(temp_dbl, history);
      scale3(Fscrit * magfs_inv, contact->fs);
    } else {
      zero3(contact->fs);
    }
  }
}

/* ----------------------------------------------------------------------
   Mindlin from pair gran/hertz/history
------------------------------------------------------------------------- */

TangentialMindlinClassic::TangentialMindlinClassic(LAMMPS *lmp) : TangentialLinearHistoryClassic(lmp)
{
  scale_area = 1; // Sets gran/hertz/history behavior
}

/* ----------------------------------------------------------------------
   Mindlin model
------------------------------------------------------------------------- */

TangentialMindlin::TangentialMindlin(LAMMPS *lmp) : TangentialModel(lmp)
{
  num_coeffs = 3;
  size_history = 3;
  mindlin_force = 0;
  mindlin_rescale = 0;
}

/* ---------------------------------------------------------------------- */

void TangentialMindlin::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu = coeffs[2];

  if (k == -1) {
    if (!contact->normal_model->material_properties)
      error->all(FLERR, "Must either specify tangential stiffness or material properties for normal model for the Mindlin tangential style");

    double Emod = contact->normal_model->Emod;
    double poiss = contact->normal_model->poiss;

    if (contact->contact_type == PAIR) {
      k = 8.0 * mix_stiffnessG(Emod, Emod, poiss, poiss);
    } else {
      k = 8.0 * mix_stiffnessG_wall(Emod, poiss);
    }
  }

  if (k < 0.0 || xt < 0.0 || mu < 0.0)
    error->all(FLERR, "Illegal Mindlin tangential model");
}

/* ---------------------------------------------------------------------- */

void TangentialMindlin::mix_coeffs(double* icoeffs, double* jcoeffs)
{
  if (icoeffs[0] == -1 || jcoeffs[0] == -1) coeffs[0] = -1;
  else coeffs[0] = mix_geom(icoeffs[0], jcoeffs[0]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

void TangentialMindlin::calculate_forces()
{
  double k_scaled, magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  damp = xt * contact->damping_model->damp_prefactor;

  double *history = & contact->history[history_index];
  double Fscrit = contact->normal_model->Fncrit * mu;

  k_scaled = k * contact->area;

  // on unloading, rescale the shear displacements/force
  if (mindlin_rescale)
    if (contact->area < history[3])
      scale3(contact->area / history[3], history);

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (contact->history_update) {
    rsht = dot3(history, contact->nx);
    if (mindlin_force) {
      frame_update = fabs(rsht) > (EPSILON * Fscrit);
    } else {
      frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
    }

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, contact->nx, temp_array);
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
      scale3(-k_scaled * contact->dt, contact->vtr, temp_array);
    } else {
      scale3(contact->dt, contact->vtr, temp_array);
    }
    add3(history, temp_array, history);

    if (mindlin_rescale) history[3] = contact->area;
  }

  // tangential forces = history + tangential velocity damping
  scale3(-damp, contact->vtr, contact->fs);

  if (!mindlin_force) {
    scale3(k_scaled, history, temp_array);
    sub3(contact->fs, temp_array, contact->fs);
  } else {
    add3(contact->fs, history, contact->fs);
  }

  // rescale frictional displacements and forces if needed
  magfs = len3(contact->fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, contact->fs, history);
      scale3(damp, contact->vtr, temp_array);
      add3(history, temp_array, history);

      if (!mindlin_force)
        scale3(-1.0 / k_scaled, history);

      scale3(Fscrit * magfs_inv, contact->fs);
    } else {
      zero3(contact->fs);
    }
  }
}

/* ----------------------------------------------------------------------
   Mindlin force model
------------------------------------------------------------------------- */

TangentialMindlinForce::TangentialMindlinForce(LAMMPS *lmp) : TangentialMindlin(lmp)
{
  num_coeffs = 3;
  size_history = 3;
  mindlin_force = 1;
  mindlin_rescale = 0;
}

/* ----------------------------------------------------------------------
   Mindlin rescale model
------------------------------------------------------------------------- */

TangentialMindlinRescale::TangentialMindlinRescale(LAMMPS *lmp) : TangentialMindlin(lmp)
{
  num_coeffs = 3;
  size_history = 4;
  mindlin_force = 0;
  mindlin_rescale = 1;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

/* ----------------------------------------------------------------------
   Mindlin rescale force model
------------------------------------------------------------------------- */

TangentialMindlinRescaleForce::TangentialMindlinRescaleForce(LAMMPS *lmp) : TangentialMindlin(lmp)
{
  num_coeffs = 3;
  size_history = 4;
  mindlin_force = 1;
  mindlin_rescale = 1;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

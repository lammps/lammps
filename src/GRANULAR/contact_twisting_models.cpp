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

#include "contact_normal_models.h"
#include "contact_tangential_models.h"
#include "contact_twisting_models.h"
#include "contact.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace Contact;
using namespace MathConst;

/* ----------------------------------------------------------------------
   Default twisting model
------------------------------------------------------------------------- */

TwistingModel::TwistingModel(LAMMPS *lmp) : SubModel(lmp) {}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

TwistingNone::TwistingNone(LAMMPS *lmp) : TwistingModel(lmp) {}

/* ----------------------------------------------------------------------
   Marshall twisting model
------------------------------------------------------------------------- */

TwistingMarshall::TwistingMarshall(LAMMPS *lmp) : TwistingModel(lmp)
{
  num_coeffs = 0;
  size_history = 3;
}

/* ---------------------------------------------------------------------- */


void TwistingMarshall::init()
{
  k_tang = contact->tangential_model->k;
  damp_tang = contact->tangential_model->damp_tangential;
  mu_tang = contact->tangential_model->mu;
}

/* ---------------------------------------------------------------------- */

void TwistingMarshall::calculate_forces()
{
  double signtwist, Mtcrit;

  // Calculate twist coefficients from tangential model & contact geometry
  // eq 32 of Marshall paper
  double k = 0.5 * k_tang * contact->area * contact->area;
  double damp = 0.5 * damp_tang * contact->area * contact->area;
  double mu = TWOTHIRDS * mu_tang * contact->area;

  if (contact->history_update) {
    contact->history[history_index] += contact->magtwist * contact->dt;
  }

  // M_t torque (eq 30)
  contact->magtortwist = -k * contact->history[history_index] - damp * contact->magtwist;
  signtwist = (contact->magtwist > 0) - (contact->magtwist < 0);
  Mtcrit = mu * contact->normal_model->Fncrit; // critical torque (eq 44)

  if (fabs(contact->magtortwist) > Mtcrit) {
    contact->history[history_index] = (Mtcrit * signtwist - damp * contact->magtwist) / k;
    contact->magtortwist = -Mtcrit * signtwist; // eq 34
  }
}

/* ----------------------------------------------------------------------
   SDS twisting model
------------------------------------------------------------------------- */

TwistingSDS::TwistingSDS(LAMMPS *lmp) : TwistingModel(lmp)
{
  num_coeffs = 3;
  size_history = 3;
}

/* ---------------------------------------------------------------------- */

void TwistingSDS::coeffs_to_local()
{
  k = coeffs[0];
  mu = coeffs[1];
  damp = coeffs[2];

  if (k < 0.0 || mu < 0.0 || damp < 0.0)
    error->all(FLERR, "Illegal SDS twisting model");
}

/* ---------------------------------------------------------------------- */

void TwistingSDS::calculate_forces()
{
  double signtwist, Mtcrit;

  if (contact->history_update) {
    contact->history[history_index] += contact->magtwist * contact->dt;
  }

  // M_t torque (eq 30)
  contact->magtortwist = -k * contact->history[history_index] - damp * contact->magtwist;
  signtwist = (contact->magtwist > 0) - (contact->magtwist < 0);
  Mtcrit = mu * contact->normal_model->Fncrit; // critical torque (eq 44)

  if (fabs(contact->magtortwist) > Mtcrit) {
    contact->history[history_index] = (Mtcrit * signtwist - damp * contact->magtwist) / k;
    contact->magtortwist = -Mtcrit * signtwist; // eq 34
  }
}

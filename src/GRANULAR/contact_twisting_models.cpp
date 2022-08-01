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
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace Contact;
using namespace MathConst;

/* ----------------------------------------------------------------------
   Marshall twisting model
------------------------------------------------------------------------- */

TwistingMarshall::TwistingMarshall()
{
  num_coeffs = 0;
  size_history = 3;
}

/* ---------------------------------------------------------------------- */

double TwistingMarshall::calculate_forces()
{
  double signtwist, Mtcrit;

  // Calculate twist coefficients from tangential model & contact geometry
  // eq 32 of Marshall paper
  double k = 0.5 * contact->tangential_model->k * contact->area * contact->area;
  double damp = 0.5 * contact->tangential_model->damp * contact->area * contact->area;
  double mu = TWOTHIRDS * contact->area * contact->tangential_model->mu;

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

TwistingSDS::TwistingSDS()
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
}

/* ---------------------------------------------------------------------- */

void TwistingSDS::mix_coeffs(TwistingModel* imodel, TwistingModel* jmodel)
{
  coeffs[0] = mix_geom(imodel->coeffs[0], jmodel->coeffs[0]);
  coeffs[1] = mix_geom(imodel->coeffs[1], jmodel->coeffs[1]);
  coeffs[2] = mix_geom(imodel->coeffs[2], jmodel->coeffs[2]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

double TwistingSDS::calculate_forces()
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

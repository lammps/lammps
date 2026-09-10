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

#include "gran_sub_mod_twisting.h"

#include "error.h"
#include "gran_sub_mod_normal.h"
#include "gran_sub_mod_tangential.h"
#include "granular_model.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace Granular_NS;

using MathConst::TWOTHIRDS;

/* ----------------------------------------------------------------------
   Default twisting model
------------------------------------------------------------------------- */

GranSubModTwisting::GranSubModTwisting(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp) {}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModTwistingNone::GranSubModTwistingNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTwisting(gm, lmp)
{
}

/* ----------------------------------------------------------------------
   Marshall twisting model
------------------------------------------------------------------------- */

GranSubModTwistingMarshall::GranSubModTwistingMarshall(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTwisting(gm, lmp)
{
  num_coeffs = 0;
  size_history = 3;
  contact_radius_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTwistingMarshall::init()
{
  k_tang = gm->tangential_model->get_k();
  mu_tang = gm->tangential_model->get_mu();
}

/* ---------------------------------------------------------------------- */

double GranSubModTwistingMarshall::calculate_forces()
{
  double signtwist, Mtcrit, magtortwist;

  double dt = gm->dt;
  double magtwist = gm->magtwist;
  double contact_radius = gm->contact_radius;
  double *history = gm->history;
  int history_update = gm->history_update;

  double Fncrit = gm->normal_model->get_fncrit();
  double tdamp = gm->tangential_model->get_damp();

  // Calculate twist coefficients from tangential model & contact geometry
  // eq 32 of Marshall paper

  double k = 0.5 * k_tang * contact_radius * contact_radius;
  double damp = 0.5 * tdamp * contact_radius * contact_radius;
  double mu = TWOTHIRDS * mu_tang * contact_radius;

  if (history_update) history[history_index] += magtwist * dt;

  // M_t torque (eq 30)
  magtortwist = -k * history[history_index] - damp * magtwist;
  signtwist = (magtwist > 0) - (magtwist < 0);
  Mtcrit = mu * Fncrit; // critical torque (eq 44)

  if (fabs(magtortwist) > Mtcrit) {
    history[history_index] = (Mtcrit * signtwist - damp * magtwist) / k;
    magtortwist = -Mtcrit * signtwist;    // eq 34
  }

  return magtortwist;
}

/* ----------------------------------------------------------------------
   SDS twisting model
------------------------------------------------------------------------- */

GranSubModTwistingSDS::GranSubModTwistingSDS(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTwisting(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
}

/* ---------------------------------------------------------------------- */

void GranSubModTwistingSDS::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];
  mu = coeffs[2];

  if (k < 0.0 || mu < 0.0 || damp < 0.0) error->all(FLERR, "Illegal SDS twisting model");
}

/* ---------------------------------------------------------------------- */

double GranSubModTwistingSDS::calculate_forces()
{
  double signtwist, Mtcrit, magtortwist;

  double magtwist = gm->magtwist;
  double dt = gm->dt;
  double *history = gm->history;
  int history_update = gm->history_update;

  double Fncrit = gm->normal_model->get_fncrit();

  if (history_update) history[history_index] += magtwist * dt;

  // M_t torque (eq 30)
  magtortwist = -k * history[history_index] - damp * magtwist;
  signtwist = (magtwist > 0) - (magtwist < 0);
  Mtcrit = mu * Fncrit;    // critical torque (eq 44)

  if (fabs(magtortwist) > Mtcrit) {
    history[history_index] = (Mtcrit * signtwist - damp * magtwist) / k;
    magtortwist = -Mtcrit * signtwist;    // eq 34
  }

  return magtortwist;
}

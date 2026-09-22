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

#include "gran_sub_mod_twisting_kernel.h"

#include "error.h"
#include "gran_sub_mod_normal.h"
#include "gran_sub_mod_tangential.h"
#include "granular_model.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace Granular_NS::GranKernel;

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

  // the twist angle is the same seen from either particle: magtwist =
  // relrot . nx is unchanged when both relrot and nx flip sign, so this
  // history must NOT be negated when it is transferred
  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = +1.0;
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
  GranTwistingParams<double> p;
  p.model = GRAN_TWISTING_NONE;    // host classes call the models directly
  p.k = 0.0;
  p.damp = 0.0;
  p.mu = 0.0;
  p.k_tang = k_tang;
  p.mu_tang = mu_tang;

  return gran_twisting_marshall(p, gm->magtwist, gm->dt, gm->contact_radius,
                                gm->normal_model->get_fncrit(),
                                gm->tangential_model->get_damp(), gm->history_update,
                                &gm->history[history_index]);
}

/* ----------------------------------------------------------------------
   SDS twisting model
------------------------------------------------------------------------- */

GranSubModTwistingSDS::GranSubModTwistingSDS(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTwisting(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;

  // see the comment in the Marshall constructor above
  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = +1.0;
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
  GranTwistingParams<double> p;
  p.model = GRAN_TWISTING_NONE;    // host classes call the models directly
  p.k = k;
  p.damp = damp;
  p.mu = mu;
  p.k_tang = 0.0;
  p.mu_tang = 0.0;

  return gran_twisting_sds(p, gm->magtwist, gm->dt, gm->normal_model->get_fncrit(),
                           gm->history_update, &gm->history[history_index]);
}

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
#include "granular_model.h"
#include "math_special.h"
#include "math_const.h"

#include <cmath>

#include <cmath>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathConst;

using MathSpecial::cube;
using MathSpecial::powint;
using MathSpecial::square;

static constexpr double TWOROOTFIVEBYSIX = 1.82574185835055380345;      // 2/sqrt(5/6)
static constexpr double ROOTTHREEBYTWO = 1.22474487139158894067;        // sqrt(3/2)

/* ----------------------------------------------------------------------
   Default damping model
------------------------------------------------------------------------- */

GranSubModDamping::GranSubModDamping(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp) {}

/* ---------------------------------------------------------------------- */

void GranSubModDamping::init()
{
  damp = gm->normal_model->get_damp();
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModDampingNone::GranSubModDampingNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModDamping(gm, lmp)
{
}

/* ---------------------------------------------------------------------- */

double GranSubModDampingNone::calculate_forces()
{
  damp_prefactor = 0.0;
  return 0.0;
}

/* ----------------------------------------------------------------------
   Velocity damping
------------------------------------------------------------------------- */

GranSubModDampingVelocity::GranSubModDampingVelocity(GranularModel *gm, LAMMPS *lmp) :
    GranSubModDamping(gm, lmp)
{
}

/* ---------------------------------------------------------------------- */

double GranSubModDampingVelocity::calculate_forces()
{
  damp_prefactor = damp;
  return -damp_prefactor * gm->vnnr;
}

/* ----------------------------------------------------------------------
   Mass velocity damping
------------------------------------------------------------------------- */

GranSubModDampingMassVelocity::GranSubModDampingMassVelocity(GranularModel *gm, LAMMPS *lmp) :
    GranSubModDamping(gm, lmp)
{
}

/* ---------------------------------------------------------------------- */

double GranSubModDampingMassVelocity::calculate_forces()
{
  damp_prefactor = damp * gm->meff;
  return -damp_prefactor * gm->vnnr;
}

/* ----------------------------------------------------------------------
   Default, viscoelastic damping
------------------------------------------------------------------------- */

GranSubModDampingViscoelastic::GranSubModDampingViscoelastic(GranularModel *gm, LAMMPS *lmp) :
    GranSubModDamping(gm, lmp)
{
  contact_radius_flag = 1;
}

/* ---------------------------------------------------------------------- */

double GranSubModDampingViscoelastic::calculate_forces()
{
  damp_prefactor = damp * gm->meff * gm->contact_radius;
  return -damp_prefactor * gm->vnnr;
}

/* ----------------------------------------------------------------------
   Tsuji damping
------------------------------------------------------------------------- */

GranSubModDampingTsuji::GranSubModDampingTsuji(GranularModel *gm, LAMMPS *lmp) :
    GranSubModDamping(gm, lmp)
{
  allow_cohesion = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModDampingTsuji::init()
{
  double tmp = gm->normal_model->get_damp();
  damp = 1.2728 - 4.2783 * tmp + 11.087 * square(tmp);
  damp += -22.348 * cube(tmp) + 27.467 * powint(tmp, 4);
  damp += -18.022 * powint(tmp, 5) + 4.8218 * powint(tmp, 6);
}

/* ---------------------------------------------------------------------- */

double GranSubModDampingTsuji::calculate_forces()
{
  // in case argument <= 0 due to precision issues
  double sqrt1;
  if (gm->delta > 0.0)
    sqrt1 = MAX(0.0, gm->meff * gm->Fnormal / gm->delta);
  else
    sqrt1 = 0.0;
  damp_prefactor = damp * sqrt(sqrt1);
  return -damp_prefactor * gm->vnnr;
}

/* ----------------------------------------------------------------------
   Coefficient of restitution damping
------------------------------------------------------------------------- */

GranSubModDampingCoeffRestitution::GranSubModDampingCoeffRestitution(GranularModel *gm, LAMMPS *lmp) :
    GranSubModDamping(gm, lmp)
{
  allow_cohesion = 0;
}

void GranSubModDampingCoeffRestitution::init()
{
  // Calculate prefactor, assume Hertzian as default
  double cor = gm->normal_model->get_damp();
  double logcor = log(cor);
  if (gm->normal_model->name == "hooke") {
    damp = -2 * logcor / sqrt(MY_PI * MY_PI + logcor * logcor);
  } else {
    damp = -ROOTTHREEBYTWO * TWOROOTFIVEBYSIX * logcor;
    damp /= sqrt(MY_PI * MY_PI + logcor * logcor);
  }
}

double GranSubModDampingCoeffRestitution::calculate_forces()
{
  // in case argument <= 0 due to precision issues
  double sqrt1;
  if (gm->delta > 0.0)
    sqrt1 = MAX(0.0, gm->meff * gm->Fnormal / gm->delta);
  else
    sqrt1 = 0.0;
  damp_prefactor = damp * sqrt(sqrt1);
  return -damp_prefactor * gm->vnnr;
}

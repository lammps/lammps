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
#include "granular_model.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathSpecial;

/* ----------------------------------------------------------------------
   Default damping model
------------------------------------------------------------------------- */

GranSubModDamping::GranSubModDamping(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp) {}

/* ---------------------------------------------------------------------- */

void GranSubModDamping::init()
{
  damp = gm->normal_model->damp;
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModDampingNone::GranSubModDampingNone(GranularModel *gm, LAMMPS *lmp) : GranSubModDamping(gm, lmp) {}

/* ---------------------------------------------------------------------- */

double GranSubModDampingNone::calculate_forces()
{
  damp_prefactor = 0.0;
  return 0.0;
}

/* ----------------------------------------------------------------------
   Velocity damping
------------------------------------------------------------------------- */

GranSubModDampingVelocity::GranSubModDampingVelocity(GranularModel *gm, LAMMPS *lmp) : GranSubModDamping(gm, lmp) {}

/* ---------------------------------------------------------------------- */

double GranSubModDampingVelocity::calculate_forces()
{
  damp_prefactor = damp;
  return -damp_prefactor * gm->vnnr;
}

/* ----------------------------------------------------------------------
   Mass velocity damping
------------------------------------------------------------------------- */

GranSubModDampingMassVelocity::GranSubModDampingMassVelocity(GranularModel *gm, LAMMPS *lmp) : GranSubModDamping(gm, lmp) {}

/* ---------------------------------------------------------------------- */

double GranSubModDampingMassVelocity::calculate_forces()
{
  damp_prefactor = damp * gm->meff;
  return -damp_prefactor * gm->vnnr;
}

/* ----------------------------------------------------------------------
   Default, viscoelastic damping
------------------------------------------------------------------------- */

GranSubModDampingViscoelastic::GranSubModDampingViscoelastic(GranularModel *gm, LAMMPS *lmp) : GranSubModDamping(gm, lmp)
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

GranSubModDampingTsuji::GranSubModDampingTsuji(GranularModel *gm, LAMMPS *lmp) : GranSubModDamping(gm, lmp)
{
  allow_cohesion = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModDampingTsuji::init()
{
  double tmp = gm->normal_model->damp;
  damp = 1.2728 - 4.2783 * tmp + 11.087 * square(tmp);
  damp += -22.348 * cube(tmp) + 27.467 * powint(tmp, 4);
  damp += -18.022 * powint(tmp, 5) + 4.8218 * powint(tmp,6);
}

/* ---------------------------------------------------------------------- */

double GranSubModDampingTsuji::calculate_forces()
{
  damp_prefactor = damp * sqrt(gm->meff * gm->Fnormal / gm->delta);
  return -damp_prefactor * gm->vnnr;
}

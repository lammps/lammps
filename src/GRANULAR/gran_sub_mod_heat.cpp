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

#include "gran_sub_mod_heat.h"

#include "error.h"
#include "granular_model.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace Granular_NS;

using MathConst::MY_PI;

/* ----------------------------------------------------------------------
   Default heat conduction
------------------------------------------------------------------------- */

GranSubModHeat::GranSubModHeat(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp) {}

/* ----------------------------------------------------------------------
   Area-based heat conduction
------------------------------------------------------------------------- */

GranSubModHeatNone::GranSubModHeatNone(GranularModel *gm, LAMMPS *lmp) : GranSubModHeat(gm, lmp) {}

/* ---------------------------------------------------------------------- */

double GranSubModHeatNone::calculate_heat()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   Radius-based heat conduction
------------------------------------------------------------------------- */

GranSubModHeatRadius::GranSubModHeatRadius(GranularModel *gm, LAMMPS *lmp) : GranSubModHeat(gm, lmp)
{
  num_coeffs = 1;
  contact_radius_flag = 1;
  conductivity = 0.0;
}

/* ---------------------------------------------------------------------- */

void GranSubModHeatRadius::coeffs_to_local()
{
  conductivity = coeffs[0];

  if (conductivity < 0.0) error->all(FLERR, "Illegal radius heat model");
}

/* ---------------------------------------------------------------------- */

double GranSubModHeatRadius::calculate_heat()
{
  return 2 * conductivity * gm->contact_radius * (gm->Tj - gm->Ti);
}


/* ----------------------------------------------------------------------
   Area-based heat conduction
------------------------------------------------------------------------- */

GranSubModHeatArea::GranSubModHeatArea(GranularModel *gm, LAMMPS *lmp) : GranSubModHeat(gm, lmp)
{
  num_coeffs = 1;
  contact_radius_flag = 1;
  heat_transfer_coeff = 0.0;
}

/* ---------------------------------------------------------------------- */

void GranSubModHeatArea::coeffs_to_local()
{
  heat_transfer_coeff = coeffs[0];

  if (heat_transfer_coeff < 0.0) error->all(FLERR, "Illegal area heat model");
}

/* ---------------------------------------------------------------------- */

double GranSubModHeatArea::calculate_heat()
{
  return heat_transfer_coeff * MY_PI * gm->contact_radius * gm->contact_radius * (gm->Tj - gm->Ti);
}

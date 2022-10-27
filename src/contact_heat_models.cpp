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

#include "contact_heat_models.h"
#include "contact.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace Contact;

/* ----------------------------------------------------------------------
   Default heat conduction
------------------------------------------------------------------------- */

HeatModel::HeatModel(LAMMPS *lmp) : SubModel(lmp) {}

/* ----------------------------------------------------------------------
   Area-based heat conduction
------------------------------------------------------------------------- */

HeatNone::HeatNone(LAMMPS *lmp) : HeatModel(lmp) {}

/* ---------------------------------------------------------------------- */

double HeatNone::calculate_heat()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   Area-based heat conduction
------------------------------------------------------------------------- */

HeatArea::HeatArea(LAMMPS *lmp) : HeatModel(lmp)
{
  num_coeffs = 1;
}

/* ---------------------------------------------------------------------- */

void HeatArea::coeffs_to_local()
{
  conductivity = coeffs[0];

  if (conductivity < 0.0) error->all(FLERR, "Illegal area heat model");
}

/* ---------------------------------------------------------------------- */

double HeatArea::calculate_heat()
{
  return conductivity * contact->area * (contact->Tj - contact->Ti);
}

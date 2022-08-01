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

using namespace LAMMPS_NS;
using namespace Contact;

/* ----------------------------------------------------------------------
   Area-based heat conduction
------------------------------------------------------------------------- */

HeatArea::HeatArea()
{
  num_coeffs = 1;
}

/* ---------------------------------------------------------------------- */

void HeatArea::coeffs_to_local()
{
  conductivity = coeffs[0];
}

/* ---------------------------------------------------------------------- */

void HeatArea::mix_coeffs(HeatModel* imodel, HeatModel* jmodel)
{
  coeffs[0] = mix_geom(imodel->coeffs[0], jmodel->coeffs[0]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

double HeatArea::calculate_heat()
{
  return conductivity * contact->area * (contact->Ti - contact->Tj);
}

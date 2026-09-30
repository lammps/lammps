/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ------------------------------------------------------
    This file is part of the BOCS package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */
#include "ldd_potential_constant.h"

#include "memory.h"
#include "utils.h"

using namespace LAMMPS_NS;

LddPotentialConstant::LddPotentialConstant(class LAMMPS *lmp) : LddPotential(lmp)
{
  n_coeffs = 1;
}

LddPotentialConstant::~LddPotentialConstant()
{
  if (allocated == 1) memory->destroy(coeffs);
  allocated = 0;
}

void LddPotentialConstant::allocate()
{
  memory->create(coeffs, n_coeffs, "ldd_potential:coeffs");
  allocated = 1;
}

void LddPotentialConstant::setup_potl(int ipt, int /*narg*/, char ** arg)
{
  if (!allocated) allocate();

  coeffs[0] = utils::numeric(FLERR, arg[ipt + 2], false, lmp);
}

double LddPotentialConstant::u(double /*rho*/)
{
  return coeffs[0];
}

double LddPotentialConstant::f(double /*rho*/)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_gyration.h"
#include "update.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGyration::ComputeGyration(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute gyration command");

  scalar_flag = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeGyration::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

double ComputeGyration::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double xcm[3];
  group->xcm(igroup,masstotal,xcm);
  scalar = group->gyration(igroup,masstotal,xcm);
  return scalar;
}

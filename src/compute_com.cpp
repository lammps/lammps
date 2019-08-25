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

#include "compute_com.h"
#include "update.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCOM::ComputeCOM(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute com command");

  vector_flag = 1;
  size_vector = 3;
  extvector = 0;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeCOM::~ComputeCOM()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeCOM::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeCOM::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (group->dynamic[igroup]) masstotal = group->mass(igroup);

  group->xcm(igroup,masstotal,vector);
}

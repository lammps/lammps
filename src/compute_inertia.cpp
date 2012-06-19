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

#include "compute_inertia.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "math_extra.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeInertia::ComputeInertia(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute inertia command");

  vector_flag = 1;
  size_vector = 3;
  extvector = 0;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

void ComputeInertia::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeInertia::compute_vector()
{
  invoked_vector = update->ntimestep;

  double xcm[3];
  group->xcm(igroup,masstotal,xcm);

  double itensor[3][3];
  group->inertia(igroup,xcm,itensor);

  double eigvect[3][3];
  int ierror;

  ierror = MathExtra::jacobi(itensor,vector,eigvect);
  if (ierror)
    error->warning(FLERR,
                   "Insufficient Jacobi rotations for principal moments of inertia");

  // sort principal moments by size
  double tmp;
  if (vector[0] < vector[1]) {
    tmp = vector[1];
    vector[1] = vector[0];
    vector[0] = tmp;
  }
  if (vector[0] < vector[2]) {
    tmp = vector[2];
    vector[2] = vector[0];
    vector[0] = tmp;
  }
  if (vector[1] < vector[2]) {
    tmp = vector[2];
    vector[2] = vector[1];
    vector[1] = tmp;
  }
  return;
}

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

#include "compute_inertia.h"

#include "error.h"
#include "group.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeInertia::ComputeInertia(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, "Illegal compute inertia command");

  vector_flag = 1;
  size_vector = 6;
  extvector = 0;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeInertia::~ComputeInertia()
{
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeInertia::init()
{
  masstotal = group->mass(igroup);
}

/* ----------------------------------------------------------------------
   compute the moment-of-inertia tensor of group of atoms around its
   center-of-mass, including the moment of inertia of finite-size particles
   the six components are ordered Ixx,Iyy,Izz,Ixy,Iyz,Ixz
------------------------------------------------------------------------- */

void ComputeInertia::compute_vector()
{
  invoked_vector = update->ntimestep;

  double xcm[3], itensor[3][3];
  if (group->dynamic[igroup]) masstotal = group->mass(igroup);
  group->xcm(igroup, masstotal, xcm);
  group->inertia(igroup, xcm, itensor);
  group->inertia_extended(igroup, itensor);

  vector[0] = itensor[0][0];
  vector[1] = itensor[1][1];
  vector[2] = itensor[2][2];
  vector[3] = itensor[0][1];
  vector[4] = itensor[1][2];
  vector[5] = itensor[0][2];
}

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

#include "mpi.h"
#include "compute_rotate_gran.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INERTIA3D 0.4          // moments of inertia for sphere and disk
#define INERTIA2D 0.5

/* ---------------------------------------------------------------------- */

ComputeRotateGran::ComputeRotateGran(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute rotate/gran command");

  if (!atom->radius_flag || !atom->rmass_flag || !atom->omega_flag)
    error->all("Compute rotate/gran requires atom attributes "
	       "radius, rmass, omega");

  scalar_flag = 1;
  extensive = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeRotateGran::init()
{
  if (domain->dimension == 3) pfactor = 0.5 * force->mvv2e * INERTIA3D;
  else pfactor = 0.5 * force->mvv2e * INERTIA2D;
}

/* ---------------------------------------------------------------------- */

double ComputeRotateGran::compute_scalar()
{
  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double erot = 0.0;

  for (int i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit)
      erot += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
	       omega[i][2]*omega[i][2]) * radius[i]*radius[i]*rmass[i];

  MPI_Allreduce(&erot,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}

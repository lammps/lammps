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
#include "compute_erotate_sphere.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INVOKED_SCALAR 1

#define INERTIA 0.4          // moment of inertia for sphere

/* ---------------------------------------------------------------------- */

ComputeERotateSphere::ComputeERotateSphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute erotate/sphere command");

  if (!atom->omega_flag) 
    error->all("Compute erotate/sphere requires atom attribute omega");

  scalar_flag = 1;
  extscalar = 1;

  inertia = new double[atom->ntypes+1];
}

/* ---------------------------------------------------------------------- */

ComputeERotateSphere::~ComputeERotateSphere()
{
  delete [] inertia;
}

/* ---------------------------------------------------------------------- */

void ComputeERotateSphere::init()
{
  pfactor = 0.5 * force->mvv2e * INERTIA;

  if (atom->mass && !atom->shape)
    error->all("Compute erotate/sphere requires atom attribute shape");
  if (!atom->mass && (!atom->radius_flag || !atom->rmass_flag))
    error->all("Compute erotate/sphere requires atom attributes radius, rmass");

  if (atom->mass) {
    double *mass = atom->mass;
    double **shape = atom->shape;
    
    for (int i = 1; i <= atom->ntypes; i++) {
      if (shape[i][0] != shape[i][1] || shape[i][0] != shape[i][2])
	error->all("Compute erotate/sphere requires spherical particle shapes");
      inertia[i] = 0.25*shape[i][0]*shape[i][0] * mass[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

double ComputeERotateSphere::compute_scalar()
{
  invoked |= INVOKED_SCALAR;

  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double erotate = 0.0;

  if (mass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	erotate += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
		    omega[i][2]*omega[i][2]) * inertia[type[i]];
  } else {
    for (int i = 0; i < nlocal; i++) 
      if (mask[i] & groupbit)
	erotate += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
		    omega[i][2]*omega[i][2]) * radius[i]*radius[i]*rmass[i];
  }

  MPI_Allreduce(&erotate,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}

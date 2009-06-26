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
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INERTIA 0.4          // moment of inertia for sphere

/* ---------------------------------------------------------------------- */

ComputeERotateSphere::ComputeERotateSphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute erotate/sphere command");

  scalar_flag = 1;
  extscalar = 1;

  // error checks

  if (!atom->omega_flag) 
    error->all("Compute erotate/sphere requires atom attribute omega");
  if (!atom->radius_flag && !atom->avec->shape_type)
    error->all("Compute erotate/sphere requires atom attribute "
	       "radius or shape");
}

/* ---------------------------------------------------------------------- */

void ComputeERotateSphere::init()
{
  int i,itype;

  // if shape used, check that all particles are spherical
  // point particles are allowed

  if (atom->radius == NULL) {
    double **shape = atom->shape;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	itype = type[i];
	if (shape[itype][0] != shape[itype][1] || 
	    shape[itype][0] != shape[itype][2])
	  error->one("Compute erotate/sphere requires "
		     "spherical particle shapes");
      }
  }

  pfactor = 0.5 * force->mvv2e * INERTIA;
}

/* ---------------------------------------------------------------------- */

double ComputeERotateSphere::compute_scalar()
{
  int i,itype;

  invoked_scalar = update->ntimestep;

  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double **shape = atom->shape;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  // sum rotational energy for each particle
  // point particles will not contribute due to radius or shape = 0

  double erotate = 0.0;

  if (radius) {
    if (rmass) {
      for (i = 0; i < nlocal; i++) 
	if (mask[i] & groupbit)
	  erotate += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
		      omega[i][2]*omega[i][2]) * radius[i]*radius[i]*rmass[i];
    } else {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	  erotate += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
		      omega[i][2]*omega[i][2]) * 
	    radius[i]*radius[i]*mass[itype];
    }

  } else {
    if (rmass) {
      for (i = 0; i < nlocal; i++) 
	if (mask[i] & groupbit) {
	  itype = type[i];
	  erotate += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
		      omega[i][2]*omega[i][2]) * 
	    shape[itype][0]*shape[itype][0]*rmass[i];
	}
    } else {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  itype = type[i];
	  erotate += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
		      omega[i][2]*omega[i][2]) *
	    shape[itype][0]*shape[itype][0]*mass[itype];
	}
    }
  }

  MPI_Allreduce(&erotate,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}

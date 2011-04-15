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
#include "compute_erotate_asphere.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeERotateAsphere::
ComputeERotateAsphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute erotate/asphere command");

  scalar_flag = 1;
  extscalar = 1;

  // error check

  if (!atom->ellipsoid_flag)
    error->all("Compute erotate/asphere requires atom style ellipsoid");
}

/* ---------------------------------------------------------------------- */

void ComputeERotateAsphere::init()
{
  // check that all particles are finite-size
  // no point particles allowed, spherical is OK

  double **shape = atom->shape;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (shape[i][0] == 0.0)
	error->one("Compute erotate/asphere requires extended particles");

  pfactor = 0.5 * force->mvv2e;
}

/* ---------------------------------------------------------------------- */

double ComputeERotateAsphere::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **quat = atom->quat;
  double **angmom = atom->angmom;
  double **shape = atom->shape;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // sum rotational energy for each particle
  // no point particles since divide by inertia

  double wbody[3],inertia[3];
  double rot[3][3];
  double erotate = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      // principal moments of inertia

      inertia[0] = rmass[i] * 
	(shape[i][1]*shape[i][1]+shape[i][2]*shape[i][2]) / 5.0;
      inertia[1] = rmass[i] * 
	(shape[i][0]*shape[i][0]+shape[i][2]*shape[i][2]) / 5.0;
      inertia[2] = rmass[i] * 
	(shape[i][0]*shape[i][0]+shape[i][1]*shape[i][1]) / 5.0;

      // wbody = angular velocity in body frame

      MathExtra::quat_to_mat(quat[i],rot);
      MathExtra::transpose_times_column3(rot,angmom[i],wbody);
      wbody[0] /= inertia[0];
      wbody[1] /= inertia[1];
      wbody[2] /= inertia[2];
      
      erotate += inertia[0]*wbody[0]*wbody[0] +
	inertia[1]*wbody[1]*wbody[1] + inertia[2]*wbody[2]*wbody[2];
    }

  MPI_Allreduce(&erotate,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}

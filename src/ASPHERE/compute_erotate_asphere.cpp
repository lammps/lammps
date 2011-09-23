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
#include "atom_vec_ellipsoid.h"
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
  if (narg != 3) error->all(FLERR,"Illegal compute erotate/asphere command");

  scalar_flag = 1;
  extscalar = 1;

  // error check

  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) 
    error->all(FLERR,"Compute erotate/asphere requires atom style ellipsoid");
}

/* ---------------------------------------------------------------------- */

void ComputeERotateAsphere::init()
{
  // check that all particles are finite-size
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
	error->one(FLERR,"Compute erotate/asphere requires extended particles");

  pfactor = 0.5 * force->mvv2e;
}

/* ---------------------------------------------------------------------- */

double ComputeERotateAsphere::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // sum rotational energy for each particle
  // no point particles since divide by inertia

  double *shape,*quat;
  double wbody[3],inertia[3];
  double rot[3][3];
  double erotate = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      shape = bonus[ellipsoid[i]].shape;
      quat = bonus[ellipsoid[i]].quat;

      // principal moments of inertia

      inertia[0] = rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]) / 5.0;
      inertia[1] = rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]) / 5.0;
      inertia[2] = rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]) / 5.0;

      // wbody = angular velocity in body frame

      MathExtra::quat_to_mat(quat,rot);
      MathExtra::transpose_matvec(rot,angmom[i],wbody);
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

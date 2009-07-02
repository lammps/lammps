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

  inertia = 
    memory->create_2d_double_array(atom->ntypes+1,3,"fix_temp_sphere:inertia");

  // error checks

  if (!atom->angmom_flag || !atom->quat_flag ||
      !atom->avec->shape_type)
    error->all("Compute erotate/asphere requires atom attributes "
	       "angmom, quat, shape");
  if (atom->radius_flag || atom->rmass_flag)
    error->all("Compute erotate/asphere cannot be used with atom attributes "
	       "diameter or rmass");
}

/* ---------------------------------------------------------------------- */

ComputeERotateAsphere::~ComputeERotateAsphere()
{
  memory->destroy_2d_double_array(inertia);
}

/* ---------------------------------------------------------------------- */

void ComputeERotateAsphere::init()
{
  // check that all particles are finite-size
  // no point particles allowed, spherical is OK

  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (shape[type[i]][0] == 0.0)
	error->one("Compute erotate/asphere requires extended particles");

  pfactor = 0.5 * force->mvv2e;
  calculate_inertia();
}

/* ---------------------------------------------------------------------- */

double ComputeERotateAsphere::compute_scalar()
{
  int i,itype;

  invoked_scalar = update->ntimestep;

  double **quat = atom->quat;
  double **angmom = atom->angmom;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  // sum rotational energy for each particle
  // no point particles since divide by inertia

  double wbody[3];
  double rot[3][3];
  double erotate = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      itype = type[i];

      // wbody = angular velocity in body frame

      MathExtra::quat_to_mat(quat[i],rot);
      MathExtra::transpose_times_column3(rot,angmom[i],wbody);
      wbody[0] /= inertia[itype][0];
      wbody[1] /= inertia[itype][1];
      wbody[2] /= inertia[itype][2];
      
      erotate += inertia[itype][0]*wbody[0]*wbody[0]+
	inertia[itype][1]*wbody[1]*wbody[1]+
	inertia[itype][2]*wbody[2]*wbody[2];
    }

  MPI_Allreduce(&erotate,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}

/* ----------------------------------------------------------------------
   principal moments of inertia for ellipsoids
------------------------------------------------------------------------- */

void ComputeERotateAsphere::calculate_inertia()
{
  double *mass = atom->mass;
  double **shape = atom->shape;

  for (int i = 1; i <= atom->ntypes; i++) {
    inertia[i][0] = mass[i] * 
      (shape[i][1]*shape[i][1]+shape[i][2]*shape[i][2]) / 5.0;
    inertia[i][1] = mass[i] * 
      (shape[i][0]*shape[i][0]+shape[i][2]*shape[i][2]) / 5.0;
    inertia[i][2] = mass[i] * 
      (shape[i][0]*shape[i][0]+shape[i][1]*shape[i][1]) / 5.0;
  }
}

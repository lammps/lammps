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
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INVOKED_SCALAR 1

/* ---------------------------------------------------------------------- */

ComputeERotateASphere::ComputeERotateASphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute erotate/asphere command");

  if (!atom->angmom_flag || !atom->quat_flag)
    error->all("Compute erotate/asphere requires atom attributes angmom, quat");

  scalar_flag = 1;
  extscalar = 1;

  inertia = 
    memory->create_2d_double_array(atom->ntypes+1,3,"fix_temp_sphere:inertia");
}

/* ---------------------------------------------------------------------- */

ComputeERotateASphere::~ComputeERotateASphere()
{
  memory->destroy_2d_double_array(inertia);
}

/* ---------------------------------------------------------------------- */

void ComputeERotateASphere::init()
{
  pfactor = 0.5 * force->mvv2e;

  if (!atom->shape)
    error->all("Compute erotate/asphere requires atom attribute shape");

  calculate_inertia();
}

/* ---------------------------------------------------------------------- */

double ComputeERotateASphere::compute_scalar()
{
  invoked |= INVOKED_SCALAR;

  double **quat = atom->quat;
  double **angmom = atom->angmom;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int itype;
  double wbody[3];
  double rot[3][3];
  double erotate = 0.0;

  for (int i = 0; i < nlocal; i++)
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

void ComputeERotateASphere::calculate_inertia()
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

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
#include "compute_rotate_dipole.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INVOKED_SCALAR 1
#define INERTIA3D 0.4          // moments of inertia for sphere and disk
#define INERTIA2D 0.5

/* ---------------------------------------------------------------------- */

ComputeRotateDipole::ComputeRotateDipole(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute rotate/dipole command");

  if (atom->dipole == NULL || !atom->omega_flag)
    error->all("Compute rotate/dipole requires atom attributes dipole, omega");

  scalar_flag = 1;
  extscalar = 1;

  inertia = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeRotateDipole::~ComputeRotateDipole()
{
  delete [] inertia;
}

/* ---------------------------------------------------------------------- */

void ComputeRotateDipole::init()
{
  delete [] inertia;
  inertia = new double[atom->ntypes+1];

  double *mass = atom->mass;
  double **shape = atom->shape;

  if (domain->dimension == 3)
    for (int i = 1; i <= atom->ntypes; i++)
      inertia[i] = INERTIA3D * mass[i] * 0.25*shape[i][0]*shape[i][0];
  else
    for (int i = 1; i <= atom->ntypes; i++)
      inertia[i] = INERTIA2D * mass[i] * 0.25*shape[i][0]*shape[i][0];
}

/* ---------------------------------------------------------------------- */

double ComputeRotateDipole::compute_scalar()
{
  invoked |= INVOKED_SCALAR;

  double *dipole = atom->dipole;
  double **omega = atom->omega;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double erot = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && dipole[type[i]] > 0.0)
      erot += inertia[type[i]] * 
	(omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
	 omega[i][2]*omega[i][2]);

  MPI_Allreduce(&erot,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= 0.5;
  return scalar;
}

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
#include "force.h"
#include "pair.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INERTIA3D 0.4          // moments of inertia for sphere and disk
#define INERTIA2D 0.5

/* ---------------------------------------------------------------------- */

ComputeRotateDipole::ComputeRotateDipole(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute rotate/dipole command");

  if (atom->check_style("dipole") == 0)
    error->all("Must use atom style dipole with compute rotate/dipole");

  scalar_flag = 1;
  extensive = 1;

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

  // insure use of dipole pair_style
  // set sigma to Pair's sigma

  Pair *pair = force->pair_match("dipole");
  if (pair == NULL)
    error->all("Pair style is incompatible with compute rotate/dipole");
  double **sigma;
  pair->extract_dipole(&sigma);

  if (force->dimension == 3)
    for (int i = 1; i <= atom->ntypes; i++)
      inertia[i] = INERTIA3D * mass[i] * 0.25*sigma[i][i]*sigma[i][i];
  else
    for (int i = 1; i <= atom->ntypes; i++)
      inertia[i] = INERTIA2D * mass[i] * 0.25*sigma[i][i]*sigma[i][i];
}

/* ---------------------------------------------------------------------- */

double ComputeRotateDipole::compute_scalar()
{
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

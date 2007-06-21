/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "compute_temp_dipole.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "error.h"

using namespace LAMMPS_NS;

// moment of inertia for a sphere

#define INERTIA 0.4

/* ---------------------------------------------------------------------- */

ComputeTempDipole::ComputeTempDipole(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute temp/dipole command");

  if (!atom->omega_flag || atom->shape == NULL)
    error->all("Compute temp/dipole requires atom attributes omega, shape");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extensive = 0;
  tempflag = 1;

  vector = new double[6];
  inertia = new double[atom->ntypes + 1];
}

/* ---------------------------------------------------------------------- */

ComputeTempDipole::~ComputeTempDipole()
{
  delete [] vector;
  delete [] inertia;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDipole::init()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  recount();

  // moment of inertia for each particle type

  double *mass = atom->mass;
  double **shape = atom->shape;

  for (int i = 1; i <= atom->ntypes; i++)
    inertia[i] = INERTIA * mass[i] * 0.25*shape[i][0]*shape[i][0];
}

/* ---------------------------------------------------------------------- */

void ComputeTempDipole::recount()
{
  double natoms = group->count(igroup);
  dof = 2.0 * force->dimension * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempDipole::compute_scalar()
{
  double **v = atom->v;
  double *mass = atom->mass;
  double **omega = atom->omega;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // rotational and translational kinetic energy

  double t = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * 
	mass[type[i]] 
	+ (omega[i][0] * omega[i][0] + omega[i][1] * omega[i][1] +
	   omega[i][2] * omega[i][2]) * inertia[type[i]];

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) recount();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDipole::compute_vector()
{
  int i;

  double **v = atom->v;
  double *mass = atom->mass;
  double **omega = atom->omega;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double rmass,imass,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  // rotational and translational kinetic energy

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      rmass = mass[type[i]];
      imass = inertia[type[i]];
      t[0] += rmass*v[i][0]*v[i][0] + imass*omega[i][0]*omega[i][0];
      t[1] += rmass*v[i][1]*v[i][1] + imass*omega[i][1]*omega[i][1];
      t[2] += rmass*v[i][2]*v[i][2] + imass*omega[i][2]*omega[i][2];
      t[3] += rmass*v[i][0]*v[i][1] + imass*omega[i][0]*omega[i][1];
      t[4] += rmass*v[i][0]*v[i][2] + imass*omega[i][0]*omega[i][2];
      t[5] += rmass*v[i][1]*v[i][2] + imass*omega[i][1]*omega[i][2];
    }

  MPI_Allreduce(&t,&vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

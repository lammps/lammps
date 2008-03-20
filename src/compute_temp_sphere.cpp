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
#include "compute_temp_sphere.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2

#define INERTIA 0.4          // moment of inertia for sphere

/* ---------------------------------------------------------------------- */

ComputeTempSphere::ComputeTempSphere(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute temp command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  vector = new double[6];
  inertia = new double[atom->ntypes+1];
}

/* ---------------------------------------------------------------------- */

ComputeTempSphere::~ComputeTempSphere()
{
  delete [] vector;
  delete [] inertia;
}

/* ---------------------------------------------------------------------- */

void ComputeTempSphere::init()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  recount();

  if (id_bias) {
    tempbias = 1;
    int i = modify->find_compute(id_bias);
    if (i < 0) error->all("Could not find compute ID for temperature bias");
    tbias = modify->compute[i];
  } else tempbias = 0;

  if (atom->mass) {
    double *mass = atom->mass;
    double **shape = atom->shape;
    
    for (int i = 1; i <= atom->ntypes; i++) {
      if (shape[i][0] != shape[i][1] || shape[i][0] != shape[i][2])
	error->all("Compute temp/sphere requires spherical particle shapes");
      inertia[i] = INERTIA * 0.25*shape[i][0]*shape[i][0] * mass[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTempSphere::recount()
{
  double natoms = group->count(igroup);
  if (domain->dimension == 3) dof = 6.0 * natoms;
  else dof = 3.0 * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempSphere::compute_scalar()
{
  invoked |= INVOKED_SCALAR;

  if (tbias) {
    if (!(tbias->invoked & INVOKED_SCALAR))
      double tmp = tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  double **v = atom->v;
  double **omega = atom->omega;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  if (mass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * 
	  mass[type[i]];
	t += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
	      omega[i][2]*omega[i][2]) * inertia[type[i]];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
	t += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
	      omega[i][2]*omega[i][2]) * INERTIA*radius[i]*radius[i]*rmass[i];
      }
  }

  if (tbias) tbias->restore_bias_all();

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) recount();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempSphere::compute_vector()
{
  int i;

  invoked |= INVOKED_VECTOR;

  if (tbias) {
    if (!(tbias->invoked & INVOKED_VECTOR)) tbias->compute_vector();
    tbias->remove_bias_all();
  }

  double **v = atom->v;
  double **omega = atom->omega;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,inertiaone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  if (mass) {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	massone = mass[type[i]];
	t[0] += massone * v[i][0]*v[i][0];
	t[1] += massone * v[i][1]*v[i][1];
	t[2] += massone * v[i][2]*v[i][2];
	t[3] += massone * v[i][0]*v[i][1];
	t[4] += massone * v[i][0]*v[i][2];
	t[5] += massone * v[i][1]*v[i][2];

	inertiaone = inertia[type[i]];
	t[0] += massone * omega[i][0]*omega[i][0];
	t[1] += massone * omega[i][1]*omega[i][1];
	t[2] += massone * omega[i][2]*omega[i][2];
	t[3] += massone * omega[i][0]*omega[i][1];
	t[4] += massone * omega[i][0]*omega[i][2];
	t[5] += massone * omega[i][1]*omega[i][2];
      }
  } else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	massone = rmass[i];
	t[0] += massone * v[i][0]*v[i][0];
	t[1] += massone * v[i][1]*v[i][1];
	t[2] += massone * v[i][2]*v[i][2];
	t[3] += massone * v[i][0]*v[i][1];
	t[4] += massone * v[i][0]*v[i][2];
	t[5] += massone * v[i][1]*v[i][2];

	inertiaone = INERTIA*radius[i]*radius[i]*rmass[i];
	t[0] += massone * omega[i][0]*omega[i][0];
	t[1] += massone * omega[i][1]*omega[i][1];
	t[2] += massone * omega[i][2]*omega[i][2];
	t[3] += massone * omega[i][0]*omega[i][1];
	t[4] += massone * omega[i][0]*omega[i][2];
	t[5] += massone * omega[i][1]*omega[i][2];
      }
  }

  if (tbias) tbias->restore_bias_all();

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempSphere::remove_bias(int i, double *v)
{
  if (tbias) tbias->remove_bias(i,v);
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempSphere::remove_bias_all()
{
  if (tbias) tbias->remove_bias_all();
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempSphere::restore_bias(double *v)
{
  if (tbias) tbias->restore_bias(v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempSphere::restore_bias_all()
{
  if (tbias) tbias->restore_bias_all();
}


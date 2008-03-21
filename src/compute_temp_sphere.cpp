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
#include "string.h"
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
  if (narg != 3 && narg != 4)
    error->all("Illegal compute temp/sphere command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  tempbias = 0;
  id_bias = NULL;
  if (narg == 4) {
    tempbias = 1;
    int n = strlen(arg[3]) + 1;
    id_bias = new char[n];
    strcpy(id_bias,arg[3]);
  }

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
  if (tempbias) {
    int i = modify->find_compute(id_bias);
    if (i < 0) error->all("Could not find compute ID for temperature bias");
    tbias = modify->compute[i];
    if (tbias->tempflag == 0)
      error->all("Bias compute does not calculate temperature");
    if (tbias->tempbias == 0)
      error->all("Bias compute does not calculate a velocity bias");
    if (tbias->igroup != igroup)
      error->all("Bias compute group does not match compute group");
    tbias->init();
    if (strcmp(tbias->style,"temp/region") == 0) tempbias = 2;
    else tempbias = 1;
  }

  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  dof_compute();

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

void ComputeTempSphere::dof_compute()
{
  double natoms = group->count(igroup);
  int nper = 6;
  if (domain->dimension == 2) nper = 3;
  dof = nper * natoms;

  if (tempbias) {
    if (tempbias == 1) dof -= tbias->dof_remove(-1) * natoms;
    else {
      int *mask = atom->mask;
      int nlocal = atom->nlocal;
      int count = 0;
      for (int i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	  if (tbias->dof_remove(i)) count++;
      int count_all;
      MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);
      dof -= nper * count_all;
    }
  }

  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempSphere::compute_scalar()
{
  invoked |= INVOKED_SCALAR;

  if (tempbias) {
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

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic || tempbias == 2) dof_compute();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempSphere::compute_vector()
{
  int i;

  invoked |= INVOKED_VECTOR;

  if (tempbias) {
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

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempSphere::remove_bias(int i, double *v)
{
  tbias->remove_bias(i,v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempSphere::restore_bias(int i, double *v)
{
  tbias->restore_bias(i,v);
}

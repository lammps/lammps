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
#include "compute_temp_region.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2

/* ---------------------------------------------------------------------- */

ComputeTempRegion::ComputeTempRegion(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all("Illegal compute temp/region command");

  iregion = domain->find_region(arg[3]);
  if (iregion == -1) error->all("Temperature region ID does not exist");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  tempbias = 1;

  maxbias = 0;
  vbiasall = NULL;
  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTempRegion::~ComputeTempRegion()
{
  memory->destroy_2d_double_array(vbiasall);
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRegion::init()
{
  dof = 0;

  tbias = NULL;
  if (id_bias) {
    int i = modify->find_compute(id_bias);
    if (i < 0) error->all("Could not find compute ID for temperature bias");
    tbias = modify->compute[i];
  }
}

/* ---------------------------------------------------------------------- */

double ComputeTempRegion::compute_scalar()
{
  invoked |= INVOKED_SCALAR;

  if (tbias) {
    if (!(tbias->invoked & INVOKED_SCALAR))
      double tmp = tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int count = 0;
  double t = 0.0;

  if (mass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit &&
	  domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2])) {
	count++;
	t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * 
	  mass[type[i]];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit &&
	  domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2])) {
	count++;
	t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
      }
  }
  
  if (tbias) tbias->restore_bias_all();

  double tarray[2],tarray_all[2];
  tarray[0] = count;
  tarray[1] = t;
  MPI_Allreduce(tarray,tarray_all,2,MPI_DOUBLE,MPI_SUM,world);
  dof = domain->dimension * tarray_all[0] - extra_dof;
  if (dof > 0) scalar = force->mvv2e * tarray_all[1] / (dof * force->boltz);
  else scalar = 0.0;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRegion::compute_vector()
{
  int i;

  invoked |= INVOKED_VECTOR;

  if (tbias) {
    if (!(tbias->invoked & INVOKED_VECTOR)) tbias->compute_vector();
    tbias->remove_bias_all();
  }

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && 
	domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2])) {
      if (mass) massone = mass[type[i]];
      else massone = rmass[i];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
    }

  if (tbias) tbias->restore_bias_all();

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRegion::remove_bias(int i, double *v)
{
  if (tbias) tbias->remove_bias(i,v);

  double *x = atom->x[i];
  if (atom->mask[i] & groupbit && 
      domain->regions[iregion]->match(x[0],x[1],x[2]))
    vbias[0] = vbias[1] = vbias[2] = 0.0;
  else {
    vbias[0] = v[0];
    vbias[1] = v[1];
    vbias[2] = v[2];
    v[0] = v[1] = v[2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRegion::remove_bias_all()
{
  if (tbias) tbias->remove_bias_all();

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (nlocal > maxbias) {
    memory->destroy_2d_double_array(vbiasall);
    maxbias = atom->nmax;
    vbiasall = memory->create_2d_double_array(maxbias,3,
					      "compute/temp:vbiasall");
  }
  
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (atom->mask[i] & groupbit && 
	  domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
	vbiasall[i][0] = vbiasall[i][1] = vbiasall[i][2] = 0.0;
      else {
	vbiasall[i][0] = v[i][0];
	vbiasall[i][1] = v[i][1];
	vbiasall[i][2] = v[i][2];
	v[i][0] = v[i][1] = v[i][2] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempRegion::restore_bias(double *v)
{
  v[0] += vbias[0];
  v[1] += vbias[1];
  v[2] += vbias[2];
  if (tbias) tbias->restore_bias(v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempRegion::restore_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] += vbiasall[i][0];
      v[i][1] += vbiasall[i][1];
      v[i][2] += vbiasall[i][2];
    }

  if (tbias) tbias->restore_bias_all();
}

/* ---------------------------------------------------------------------- */

double ComputeTempRegion::memory_usage()
{
  double bytes = maxbias * sizeof(double);
  return bytes;
}

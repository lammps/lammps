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

#include "compute_temp_region.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempRegion::ComputeTempRegion(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idregion(NULL)
{
  if (narg != 4) error->all(FLERR,"Illegal compute temp/region command");

  iregion = domain->find_region(arg[3]);
  if (iregion == -1)
    error->all(FLERR,"Region ID for compute temp/region does not exist");
  int n = strlen(arg[3]) + 1;
  idregion = new char[n];
  strcpy(idregion,arg[3]);

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
  delete [] idregion;
  memory->destroy(vbiasall);
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRegion::init()
{
  // set index and check validity of region

  iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all(FLERR,"Region ID for compute temp/region does not exist");
}

/* ---------------------------------------------------------------------- */

void ComputeTempRegion::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof = 0.0;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRegion::dof_remove_pre()
{
  domain->regions[iregion]->prematch();
}

/* ---------------------------------------------------------------------- */

int ComputeTempRegion::dof_remove(int i)
{
  double *x = atom->x[i];
  if (domain->regions[iregion]->match(x[0],x[1],x[2])) return 0;
  return 1;
}

/* ---------------------------------------------------------------------- */

double ComputeTempRegion::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  Region *region = domain->regions[iregion];
  region->prematch();

  int count = 0;
  double t = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        count++;
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        count++;
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          mass[type[i]];
      }
  }

  double tarray[2],tarray_all[2];
  tarray[0] = count;
  tarray[1] = t;
  MPI_Allreduce(tarray,tarray_all,2,MPI_DOUBLE,MPI_SUM,world);
  dof = domain->dimension * tarray_all[0] - extra_dof;
  if (dof < 0.0 && tarray_all[0] > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  if (dof > 0) scalar = force->mvv2e * tarray_all[1] / (dof * force->boltz);
  else scalar = 0.0;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRegion::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  Region *region = domain->regions[iregion];
  region->prematch();

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRegion::remove_bias(int i, double *v)
{
  double *x = atom->x[i];
  if (domain->regions[iregion]->match(x[0],x[1],x[2]))
    vbias[0] = vbias[1] = vbias[2] = 0.0;
  else {
    vbias[0] = v[0];
    vbias[1] = v[1];
    vbias[2] = v[2];
    v[0] = v[1] = v[2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRegion::remove_bias_thr(int i, double *v, double *b)
{
  double *x = atom->x[i];
  if (domain->regions[iregion]->match(x[0],x[1],x[2]))
    b[0] = b[1] = b[2] = 0.0;
  else {
    b[0] = v[0];
    b[1] = v[1];
    b[2] = v[2];
    v[0] = v[1] = v[2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRegion::remove_bias_all()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    memory->destroy(vbiasall);
    maxbias = atom->nmax;
    memory->create(vbiasall,maxbias,3,"temp/region:vbiasall");
  }

  Region *region = domain->regions[iregion];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (region->match(x[i][0],x[i][1],x[i][2]))
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

void ComputeTempRegion::restore_bias(int /*i*/, double *v)
{
  v[0] += vbias[0];
  v[1] += vbias[1];
  v[2] += vbias[2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias_thr()
   assume remove_bias_thr() was previously called with the same buffer b
------------------------------------------------------------------------- */

void ComputeTempRegion::restore_bias_thr(int /*i*/, double *v, double *b)
{
  v[0] += b[0];
  v[1] += b[1];
  v[2] += b[2];
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
}

/* ---------------------------------------------------------------------- */

double ComputeTempRegion::memory_usage()
{
  double bytes = 3*maxbias * sizeof(double);
  return bytes;
}

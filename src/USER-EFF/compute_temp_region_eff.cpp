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

/* ----------------------------------------------------------------------
   Contributing author: Andres Jaramillo-Botero (Caltech)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "compute_temp_region_eff.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempRegionEff::ComputeTempRegionEff(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (!atom->electron_flag)
    error->all(FLERR,"Compute temp/region/eff requires atom style electron");

  if (narg != 4) error->all(FLERR,"Illegal compute temp/region/eff command");

  iregion = domain->find_region(arg[3]);
  if (iregion == -1)
    error->all(FLERR,"Region ID for compute temp/region/eff does not exist");
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

ComputeTempRegionEff::~ComputeTempRegionEff()
{
  delete [] idregion;
  memory->destroy(vbiasall);
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRegionEff::init()
{
  // set index and check validity of region

  iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all(FLERR,"Region ID for compute temp/region/eff does not exist");
}

/* ---------------------------------------------------------------------- */

void ComputeTempRegionEff::setup()
{
  dof = 0.0;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRegionEff::dof_remove_pre()
{
  domain->regions[iregion]->prematch();
}

/* ---------------------------------------------------------------------- */

int ComputeTempRegionEff::dof_remove(int i)
{
  double *x = atom->x[i];
  if (domain->regions[iregion]->match(x[0],x[1],x[2])) return 0;
  return 1;
}

/* ---------------------------------------------------------------------- */

double ComputeTempRegionEff::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *ervel = atom->ervel;
  double *mass = atom->mass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double mefactor = domain->dimension/4.0;

  Region *region = domain->regions[iregion];
  region->prematch();

  int count = 0;
  int ecount = 0;
  double t = 0.0;

  if (mass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        count++;
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          mass[type[i]];
        if (fabs(spin[i])==1) {
          t += mefactor*mass[type[i]]*ervel[i]*ervel[i];
          ecount++;
        }
      }
  }

  double tarray[2],tarray_all[2];
  // Assume 3/2 k T per nucleus
  tarray[0] = count-ecount;
  tarray[1] = t;
  MPI_Allreduce(tarray,tarray_all,2,MPI_DOUBLE,MPI_SUM,world);
  dof = domain->dimension * tarray_all[0] - extra_dof;

  int one = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      if (fabs(spin[i])==1) one++;
    }

  if (dof > 0) scalar = force->mvv2e * tarray_all[1] / (dof * force->boltz);
  else scalar = 0.0;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRegionEff::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *ervel = atom->ervel;
  double *mass = atom->mass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double mefactor = domain->dimension/4.0;

  Region *region = domain->regions[iregion];
  region->prematch();

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      massone = mass[type[i]];

      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];

      if (fabs(spin[i])==1) {
        t[0] += mefactor * massone * ervel[i]*ervel[i];
        t[1] += mefactor * massone * ervel[i]*ervel[i];
        t[2] += mefactor * massone * ervel[i]*ervel[i];
      }
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
   NOTE: the following commands do not bias the radial electron velocities
------------------------------------------------------------------------- */

void ComputeTempRegionEff::remove_bias(int i, double *v)
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
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempRegionEff::remove_bias_all()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (nlocal > maxbias) {
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

void ComputeTempRegionEff::restore_bias(int i, double *v)
{
  v[0] += vbias[0];
  v[1] += vbias[1];
  v[2] += vbias[2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempRegionEff::restore_bias_all()
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

double ComputeTempRegionEff::memory_usage()
{
  double bytes = maxbias * sizeof(double);
  return bytes;
}

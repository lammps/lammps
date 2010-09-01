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
#include "string.h"
#include "compute_temp_region_eff.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempRegionEff::
ComputeTempRegionEff(LAMMPS *lmp, int narg, char **arg) : 
  ComputeTempRegion(lmp, narg, arg)
{
  // error check

  if (!atom->spin_flag || !atom->ervel_flag) 
    error->all("Compute temp/region/eff requires atom attributes spin, ervel");
}

/* ---------------------------------------------------------------------- */

double ComputeTempRegionEff::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *ervel = atom->ervel;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  Region *region = domain->regions[iregion];
  int count = 0;
  double t = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
	count++;
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2])*rmass[i];
        if (spin[i]) t += 0.75*rmass[i]*ervel[i]*ervel[i];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
	count++;
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * 
	  mass[type[i]];
        if (spin[i]) t += 0.75*mass[type[i]]*ervel[i]*ervel[i];
      }
  }

  double tarray[2],tarray_all[2];
  tarray[0] = count;
  tarray[1] = t;
  MPI_Allreduce(tarray,tarray_all,2,MPI_DOUBLE,MPI_SUM,world);
  dof = domain->dimension * tarray_all[0] - extra_dof;
 
  int one = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      if (spin[i]) one++;
    }
  int nelectrons_region;
  MPI_Allreduce(&one,&nelectrons_region,1,MPI_INT,MPI_SUM,world);

  // average over nuclear dof only

  dof -= domain->dimension * nelectrons_region ;

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
  double *rmass = atom->rmass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  Region *region = domain->regions[iregion];
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
		
      if (spin[i]) {
        t[0] += 0.75 * massone * ervel[i]*ervel[i];
        t[1] += 0.75 * massone * ervel[i]*ervel[i];
        t[2] += 0.75 * massone * ervel[i]*ervel[i];
      }			
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

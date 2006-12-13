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
#include "string.h"
#include "temp_region.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

TempRegion::TempRegion(int narg, char **arg) : Temperature(narg, arg)
{
  options(narg-4,&arg[4]);

  for (iregion = 0; iregion < domain->nregion; iregion++)
    if (strcmp(arg[3],domain->regions[iregion]->id) == 0) break;
  if (iregion == domain->nregion)
    error->all("Temperature region ID does not exist");
}

/* ---------------------------------------------------------------------- */

void TempRegion::init()
{
  dof = 0;
}

/* ---------------------------------------------------------------------- */

double TempRegion::compute()
{
  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int count = 0;
  double t = 0.0;

  if (atom->mass_require) {
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
  
  double tarray[2],tarray_all[2];
  tarray[0] = count;
  tarray[1] = t;
  MPI_Allreduce(tarray,tarray_all,2,MPI_DOUBLE,MPI_SUM,world);
  dof = force->dimension * tarray_all[0] - extra_dof;
  if (dof > 0) t_total = force->mvv2e * tarray_all[1] / (dof * force->boltz);
  else t_total = 0.0;
  return t_total;
}

/* ---------------------------------------------------------------------- */

void TempRegion::tensor()
{
  int i;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int mass_require = atom->mass_require;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && 
	domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2])) {
      if (mass_require) massone = mass[type[i]];
      else massone = rmass[i];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
    }

  MPI_Allreduce(&t,&ke_tensor,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) ke_tensor[i] *= force->mvv2e;
}

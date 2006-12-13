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
#include "stdlib.h"
#include "temp_partial.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

TempPartial::TempPartial(int narg, char **arg) : Temperature(narg, arg)
{
  options(narg-6,&arg[6]);

  xflag = atoi(arg[3]);
  yflag = atoi(arg[4]);
  zflag = atoi(arg[5]);
}

/* ---------------------------------------------------------------------- */

void TempPartial::init()
{
  count_fix();
  recount();
}

/* ---------------------------------------------------------------------- */

void TempPartial::recount()
{
  double natoms = group->count(igroup);
  dof = (xflag+yflag+zflag) * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double TempPartial::compute()
{
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  if (atom->mass_require) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	t += (xflag*v[i][0]*v[i][0] + yflag*v[i][1]*v[i][1] + 
	      zflag*v[i][2]*v[i][2]) * mass[type[i]];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	t += (xflag*v[i][0]*v[i][0] + yflag*v[i][1]*v[i][1] + 
	      zflag*v[i][2]*v[i][2]) * rmass[i];
  }
    
  MPI_Allreduce(&t,&t_total,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) recount();
  t_total *= tfactor;
  return t_total;
}

/* ---------------------------------------------------------------------- */

void TempPartial::tensor()
{
  int i;

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
    if (mask[i] & groupbit) {
      if (mass_require) massone = mass[type[i]];
      else massone = rmass[i];
      t[0] += massone * xflag*v[i][0]*v[i][0];
      t[1] += massone * yflag*v[i][1]*v[i][1];
      t[2] += massone * zflag*v[i][2]*v[i][2];
      t[3] += massone * xflag*yflag*v[i][0]*v[i][1];
      t[4] += massone * xflag*zflag*v[i][0]*v[i][2];
      t[5] += massone * yflag*zflag*v[i][1]*v[i][2];
    }

  MPI_Allreduce(&t,&ke_tensor,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) ke_tensor[i] *= force->mvv2e;
}

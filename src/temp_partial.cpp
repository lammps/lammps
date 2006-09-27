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
  xflag = atoi(arg[3]);
  yflag = atoi(arg[4]);
  zflag = atoi(arg[5]);
}

/* ---------------------------------------------------------------------- */

void TempPartial::init()
{
  count_atoms();
  count_fix();
  dof = (xflag+yflag+zflag) * ncount;
  dof -= extra_dof + fix_dof;
  if (ncount > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double TempPartial::compute()
{
  double **v = atom->v;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      t += (xflag*v[i][0]*v[i][0] + yflag*v[i][1]*v[i][1] + 
	    zflag*v[i][2]*v[i][2]) * mass[type[i]];

  MPI_Allreduce(&t,&t_total,1,MPI_DOUBLE,MPI_SUM,world);
  t_total *= tfactor;
  return t_total;
}

/* ---------------------------------------------------------------------- */

void TempPartial::tensor()
{
  int i;

  double **v = atom->v;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double rmass,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      rmass = mass[type[i]];
      t[0] += rmass * xflag*v[i][0]*v[i][0];
      t[1] += rmass * yflag*v[i][1]*v[i][1];
      t[2] += rmass * zflag*v[i][2]*v[i][2];
      t[3] += rmass * xflag*yflag*v[i][0]*v[i][1];
      t[4] += rmass * xflag*zflag*v[i][0]*v[i][2];
      t[5] += rmass * yflag*zflag*v[i][1]*v[i][2];
    }

  MPI_Allreduce(&t,&ke_tensor,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) ke_tensor[i] *= force->mvv2e;
}

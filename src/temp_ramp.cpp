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
#include "string.h"
#include "temp_ramp.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "error.h"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

TempRamp::TempRamp(int narg, char **arg) : Temperature(narg, arg)
{
  options(narg-9,&arg[9]);

  if (strcmp(arg[3],"vx") == 0) v_dim = 0;
  else if (strcmp(arg[3],"vy") == 0) v_dim = 1;
  else if (strcmp(arg[3],"vz") == 0) v_dim = 2;
  else error->all("Illegal temperature ramp command");

  if (v_dim == 0) {
    v_lo = xscale*atof(arg[4]);
    v_hi = xscale*atof(arg[5]);
  } else if (v_dim == 1) {
    v_lo = yscale*atof(arg[4]);
    v_hi = yscale*atof(arg[5]);
  } else if (v_dim == 0) {
    v_lo = zscale*atof(arg[4]);
    v_hi = zscale*atof(arg[5]);
  }

  if (strcmp(arg[6],"x") == 0) coord_dim = 0;
  else if (strcmp(arg[6],"y") == 0) coord_dim = 1;
  else if (strcmp(arg[6],"z") == 0) coord_dim = 2;
  else error->all("Illegal temperature ramp command");

  if (coord_dim == 0) {
    coord_lo = xscale*atof(arg[7]);
    coord_hi = xscale*atof(arg[8]);
  } else if (coord_dim == 1) {
    coord_lo = yscale*atof(arg[7]);
    coord_hi = yscale*atof(arg[8]);
  } else if (coord_dim == 2) {
    coord_lo = zscale*atof(arg[7]);
    coord_hi = zscale*atof(arg[8]);
  }
}

/* ---------------------------------------------------------------------- */

void TempRamp::init()
{
  count_fix();
  recount();
}

/* ---------------------------------------------------------------------- */

void TempRamp::recount()
{
  double natoms = group->count(igroup);
  dof = force->dimension * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double TempRamp::compute()
{
  double fraction,vramp,vtmp[3];

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int mass_require = atom->mass_require;

  double t = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi - v_lo);
      vtmp[0] = v[i][0];
      vtmp[1] = v[i][1];
      vtmp[2] = v[i][2];
      vtmp[v_dim] -= vramp;
      if (mass_require)
	t += (vtmp[0]*vtmp[0] + vtmp[1]*vtmp[1] + vtmp[2]*vtmp[2]) * 
	  mass[type[i]];
      else
	t += (vtmp[0]*vtmp[0] + vtmp[1]*vtmp[1] + vtmp[2]*vtmp[2]) * rmass[i];
    }

  MPI_Allreduce(&t,&t_total,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) recount();
  t_total *= tfactor;
  return t_total;
}

/* ---------------------------------------------------------------------- */

void TempRamp::tensor()
{
  int i;
  double fraction,vramp,vtmp[3];

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
    if (mask[i] & groupbit) {
      fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi - v_lo);
      vtmp[0] = v[i][0];
      vtmp[1] = v[i][1];
      vtmp[2] = v[i][2];
      vtmp[v_dim] -= vramp;

      if (mass_require) massone = mass[type[i]];
      else massone = rmass[i];
      t[0] += massone * vtmp[0]*vtmp[0];
      t[1] += massone * vtmp[1]*vtmp[1];
      t[2] += massone * vtmp[2]*vtmp[2];
      t[3] += massone * vtmp[0]*vtmp[1];
      t[4] += massone * vtmp[0]*vtmp[2];
      t[5] += massone * vtmp[1]*vtmp[2];
    }

  MPI_Allreduce(&t,&ke_tensor,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) ke_tensor[i] *= force->mvv2e;
}

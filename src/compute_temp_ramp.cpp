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
#include "stdlib.h"
#include "string.h"
#include "compute_temp_ramp.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "lattice.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2

/* ---------------------------------------------------------------------- */

ComputeTempRamp::ComputeTempRamp(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg < 9) error->all("Illegal compute temp command");

  // parse optional args

  scaleflag = 1;

  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal compute temp/ramp command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal compute temp/ramp command");
      iarg += 2;
    } else error->all("Illegal compute temp/ramp command");
  }

  // setup scaling

  if (scaleflag && domain->lattice == NULL)
    error->all("Use of compute temp/ramp with undefined lattice");

  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // read standard args and apply scaling

  if (strcmp(arg[3],"vx") == 0) v_dim = 0;
  else if (strcmp(arg[3],"vy") == 0) v_dim = 1;
  else if (strcmp(arg[3],"vz") == 0) v_dim = 2;
  else error->all("Illegal compute temp/ramp command");

  if (v_dim == 0) {
    v_lo = xscale*atof(arg[4]);
    v_hi = xscale*atof(arg[5]);
  } else if (v_dim == 1) {
    v_lo = yscale*atof(arg[4]);
    v_hi = yscale*atof(arg[5]);
  } else if (v_dim == 2) {
    v_lo = zscale*atof(arg[4]);
    v_hi = zscale*atof(arg[5]);
  }

  if (strcmp(arg[6],"x") == 0) coord_dim = 0;
  else if (strcmp(arg[6],"y") == 0) coord_dim = 1;
  else if (strcmp(arg[6],"z") == 0) coord_dim = 2;
  else error->all("Illegal compute temp/ramp command");

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

  // settings

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTempRamp::~ComputeTempRamp()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRamp::init()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  recount();
}

/* ---------------------------------------------------------------------- */

void ComputeTempRamp::recount()
{
  double natoms = group->count(igroup);
  dof = domain->dimension * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempRamp::compute_scalar()
{
  double fraction,vramp,vthermal[3];

  invoked |= INVOKED_SCALAR;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi - v_lo);
      vthermal[0] = v[i][0];
      vthermal[1] = v[i][1];
      vthermal[2] = v[i][2];
      vthermal[v_dim] -= vramp;
      if (mass)
	t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
	      vthermal[2]*vthermal[2]) * mass[type[i]];
      else
	t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
	      vthermal[2]*vthermal[2]) * rmass[i];
    }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) recount();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempRamp::compute_vector()
{
  int i;
  double fraction,vramp,vthermal[3];

  invoked |= INVOKED_VECTOR;

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
    if (mask[i] & groupbit) {
      fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi - v_lo);
      vthermal[0] = v[i][0];
      vthermal[1] = v[i][1];
      vthermal[2] = v[i][2];
      vthermal[v_dim] -= vramp;

      if (mass) massone = mass[type[i]];
      else massone = rmass[i];
      t[0] += massone * vthermal[0]*vthermal[0];
      t[1] += massone * vthermal[1]*vthermal[1];
      t[2] += massone * vthermal[2]*vthermal[2];
      t[3] += massone * vthermal[0]*vthermal[1];
      t[4] += massone * vthermal[0]*vthermal[2];
      t[5] += massone * vthermal[1]*vthermal[2];
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

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

/* ----------------------------------------------------------------------
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "compute_temp_deform.h"
#include "domain.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "group.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NO_REMAP,X_REMAP,V_REMAP};                   // same as fix_deform.cpp

/* ---------------------------------------------------------------------- */

ComputeTempDeform::ComputeTempDeform(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute temp/deform command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extensive = 0;
  tempflag = 1;

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTempDeform::~ComputeTempDeform()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeform::init()
{
  int i;

  fix_dof = 0;
  for (i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  recount();

  // check fix deform remap settings

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      if (((FixDeform *) modify->fix[i])->remapflag == X_REMAP && 
	  comm->me == 0)
	error->warning("Using compute temp/deform with inconsistent fix deform remap option");
      break;
    }
  if (i == modify->nfix && comm->me == 0)
    error->warning("Using compute temp/deform with no fix deform defined");
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeform::recount()
{
  double natoms = group->count(igroup);
  dof = domain->dimension * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempDeform::compute_scalar()
{
  double lamda[3],vstream[3],vthermal[3];

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // lamda = 0-1 triclinic lamda coords
  // vstream = streaming velocity = Hrate*lamda + Hratelo
  // vthermal = thermal velocity = v - vstream

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  double t = 0.0;

  if (mass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	domain->x2lamda(x[i],lamda);
	vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] + 
	  h_rate[4]*lamda[2] + h_ratelo[0];
	vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
	vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
	vthermal[0] = v[i][0] - vstream[0];
	vthermal[1] = v[i][1] - vstream[1];
	vthermal[2] = v[i][2] - vstream[2];
	t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
	      vthermal[2]*vthermal[2]) * mass[type[i]];
      }
  }
  else
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	domain->x2lamda(x[i],lamda);
	vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] + 
	  h_rate[4]*lamda[2] + h_ratelo[0];
	vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
	vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
	vthermal[0] = v[i][0] - vstream[0];
	vthermal[1] = v[i][1] - vstream[1];
	vthermal[2] = v[i][2] - vstream[2];
	t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
	      vthermal[2]*vthermal[2]) * rmass[i];
      }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeform::compute_vector()
{
  double lamda[3],vstream[3],vthermal[3];

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  double massone,t[6];
  for (int i = 0; i < 6; i++) t[i] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] + 
	h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      vthermal[0] = v[i][0] - vstream[0];
      vthermal[1] = v[i][1] - vstream[1];
      vthermal[2] = v[i][2] - vstream[2];

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
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

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
#include "compute_temp_com.h"
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

ComputeTempCOM::ComputeTempCOM(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute temp command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTempCOM::~ComputeTempCOM()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempCOM::init()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  recount();
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeTempCOM::recount()
{
  double natoms = group->count(igroup);
  dof = domain->dimension * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempCOM::compute_scalar()
{
  double vcm[3],vthermal[3];

  invoked |= INVOKED_SCALAR;

  if (dynamic) masstotal = group->mass(igroup);
  group->vcm(igroup,masstotal,vcm);

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      vthermal[0] = v[i][0] - vcm[0];
      vthermal[1] = v[i][1] - vcm[1];
      vthermal[2] = v[i][2] - vcm[2];
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

void ComputeTempCOM::compute_vector()
{
  int i;
  double vcm[3],vthermal[3];

  invoked |= INVOKED_VECTOR;

  if (dynamic) masstotal = group->mass(igroup);
  group->vcm(igroup,masstotal,vcm);

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
      vthermal[0] = v[i][0] - vcm[0];
      vthermal[1] = v[i][1] - vcm[1];
      vthermal[2] = v[i][2] - vcm[2];

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

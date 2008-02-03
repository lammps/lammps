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
#include "compute_temp_com.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2

/* ---------------------------------------------------------------------- */

ComputeTempCoM::ComputeTempCoM(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute temp/com command");

  scalar_flag = vector_flag = 1;
  size_vector = 2;
  extscalar = 0;
  extvector = 0;
  tempflag = 1;

  vector = new double[2];

  // set comm size needed by this Compute

  comm_forward = 3;
}

/* ---------------------------------------------------------------------- */

ComputeTempCoM::~ComputeTempCoM()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempCoM::init()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  recount();
}

/* ---------------------------------------------------------------------- */

void ComputeTempCoM::recount()
{
  double natoms = group->count(igroup);
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;

  int atom1,m,i,nsum,nsumall;

  nsum = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      nsum += num_bond[i];
    }

  MPI_Allreduce(&nsum,&nsumall,1,MPI_INT,MPI_SUM,world);

  dof = 3 * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;

  dof1 = 3 * nsumall;
  dof1 -= extra_dof + fix_dof;
  if (dof1 > 0) tfactor1 = force->mvv2e / (dof1 * force->boltz);
  else tfactor1 = 0.0;

}

/* ---------------------------------------------------------------------- */

double ComputeTempCoM::compute_scalar()
{
  invoked |= INVOKED_SCALAR;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int *tag = atom->tag;

  double massone,t;
  double psum[3],msum;
  int atom1,m,i;

  // communicate ghost particle velocities

  comm->comm_compute(this);

  t = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (mass) massone = mass[type[i]];
      else massone = rmass[i];
      psum[0] = v[i][0] * massone;
      psum[1] = v[i][1] * massone;
      psum[2] = v[i][2] * massone;
      msum = massone;
      for (m = 0; m < num_bond[i]; m++) {
	atom1 = atom->map(bond_atom[i][m]);
	if (atom1 == -1) {
	  char str[128];
	  sprintf(str,"Bond atoms %d %d missing",
		  tag[i],bond_atom[i][m]);
	  error->one(str);
	}
	if (mass) massone = mass[type[atom1]];
	else massone = rmass[atom1];
	psum[0] += v[atom1][0] * massone;
	psum[1] += v[atom1][1] * massone;
	psum[2] += v[atom1][2] * massone;
	msum += massone;
      }
      t += (psum[0]*psum[0] + psum[1]*psum[1] + 
	    psum[2]*psum[2]) / msum;

    }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) recount();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempCoM::compute_vector()
{
  int i;

  invoked |= INVOKED_VECTOR;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int *tag = atom->tag;

  double massone,t[2];
  double psum[3],msum,vav[3];
  int atom1,m;

  // communicate ghost particle velocities

  comm->comm_compute(this);

  t[0] = 0.0;
  t[1] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (mass) massone = mass[type[i]];
      else massone = rmass[i];
      psum[0] = v[i][0] * massone;
      psum[1] = v[i][1] * massone;
      psum[2] = v[i][2] * massone;
      msum = massone;
      for (m = 0; m < num_bond[i]; m++) {
	atom1 = atom->map(bond_atom[i][m]);
	if (atom1 == -1) {
	  char str[128];
	  sprintf(str,"Bond atoms %d %d missing",
		  tag[i],bond_atom[i][m]);
	  error->one(str);
	}
	if (mass) massone = mass[type[atom1]];
	else massone = rmass[atom1];
	psum[0] += v[atom1][0] * massone;
	psum[1] += v[atom1][1] * massone;
	psum[2] += v[atom1][2] * massone;
	msum += massone;
      }
      t[0] += (psum[0]*psum[0] + psum[1]*psum[1] + 
	       psum[2]*psum[2]) / msum;

      vav[0] = psum[0]/msum;
      vav[1] = psum[1]/msum;
      vav[2] = psum[2]/msum;

      if (mass) massone = mass[type[i]];
      else massone = rmass[i];
      t[1] += ( 
	       (v[i][0]-vav[0])*(v[i][0]-vav[0]) +
	       (v[i][1]-vav[1])*(v[i][1]-vav[1]) +
	       (v[i][2]-vav[2])*(v[i][2]-vav[2]) 
		) * massone;
      for (m = 0; m < num_bond[i]; m++) {
	atom1 = atom->map(bond_atom[i][m]);
	if (atom1 == -1) {
	  char str[128];
	  sprintf(str,"Bond atoms %d %d missing",
		  tag[i],bond_atom[i][m]);
	  error->one(str);
	}
	if (mass) massone = mass[type[atom1]];
	else massone = rmass[atom1];
	t[1] += ( 
		 (v[atom1][0]-vav[0])*(v[atom1][0]-vav[0]) +
		 (v[atom1][1]-vav[1])*(v[atom1][1]-vav[1]) +
		 (v[atom1][2]-vav[2])*(v[atom1][2]-vav[2]) 
		  ) * massone;
      }

    }

  MPI_Allreduce(t,vector,2,MPI_DOUBLE,MPI_SUM,world);
  vector[0] *= tfactor;
  vector[1] *= tfactor1;
}

/* ---------------------------------------------------------------------- */

int ComputeTempCoM::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  double **v = atom->v;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = v[j][0];
    buf[m++] = v[j][1];
    buf[m++] = v[j][2];
  }
  return 3;
}

/* ---------------------------------------------------------------------- */

void ComputeTempCoM::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  double **v = atom->v;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

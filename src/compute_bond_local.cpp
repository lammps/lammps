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

#include "math.h"
#include "string.h"
#include "compute_bond_local.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "bond.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{DISTANCE,ENERGY};

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputeBondLocal::ComputeBondLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal compute bond/local command");

  if (atom->avec->bonds_allow == 0)
    error->all("Compute bond/local used when bonds are not allowed");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  which = new int[nvalues];
  pack_choice = new FnPtrPack[nvalues];

  dflag = eflag = 0;

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;

    if (strcmp(arg[iarg],"distance") == 0) {
      dflag = 1;
      which[i] = DISTANCE;
      pack_choice[i] = &ComputeBondLocal::pack_distance;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      eflag = 1;
      which[i] = ENERGY;
      pack_choice[i] = &ComputeBondLocal::pack_energy;
    } else error->all("Invalid keyword in compute bond/local command");
  }

  nmax = 0;
  distance = energy = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBondLocal::~ComputeBondLocal()
{
  delete [] which;
  delete [] pack_choice;
  memory->sfree(distance);
  memory->sfree(energy);
  memory->destroy_2d_double_array(array);
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::init()
{
  if (force->bond == NULL) 
    error->all("No bond style is defined for compute bond/local");

  // do initial memory allocation so that memory_usage() is correct

  ncount = compute_bonds(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute bond info

  ncount = compute_bonds(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_bonds(1);

  // fill array with distance/energy values

  if (nvalues > 1) {
    if (array) buf = array[0];
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   count bonds and compute bond info on this proc
   only count bond once if newton_bond is off
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if bond is deleted (type = 0), do not count
   if bond is turned off (type < 0), still count
   if flag is set, compute requested info about bond
------------------------------------------------------------------------- */

int ComputeBondLocal::compute_bonds(int flag)
{
  int i,atom1,atom2;
  double delx,dely,delz,rsq;

  double **x = atom->x;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  Bond *bond = force->bond;

  int m = 0;
  for (atom1 = 0; atom1 < nlocal; atom1++) {
    if (!(mask[atom1] & groupbit)) continue;
    for (i = 0; i < num_bond[atom1]; i++) {
      atom2 = atom->map(bond_atom[atom1][i]);
      if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;
      if (newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
      if (bond_type[atom1][i] == 0) continue;

      if (flag) {
	delx = x[atom1][0] - x[atom2][0];
	dely = x[atom1][1] - x[atom2][1];
	delz = x[atom1][2] - x[atom2][2];
	domain->minimum_image(delx,dely,delz);
	rsq = delx*delx + dely*dely + delz*delz;
	if (dflag) distance[m] = sqrt(rsq);
	if (eflag)
	  energy[m] = bond->single(bond_type[atom1][i],rsq,atom1,atom2);
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::pack_distance(int n)
{
  for (int m = 0; m < ncount; m++) {
    buf[n] = distance[m];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::pack_energy(int n)
{
  for (int m = 0; m < ncount; m++) {
    buf[n] = energy[m];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (dflag) {
    memory->sfree(distance);
    distance = (double *) memory->smalloc(nmax*sizeof(double),
					  "bond/local:distance");
  }
  if (eflag) {
    memory->sfree(energy);
    energy = (double *) memory->smalloc(nmax*sizeof(double),
					"bond/local:energy");
  }

  if (nvalues == 1) {
    if (dflag) vector_local = distance;
    if (eflag) vector_local = energy;
  } else {
    memory->destroy_2d_double_array(array);
    array = memory->create_2d_double_array(nmax,nvalues,
					   "bond/local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeBondLocal::memory_usage()
{
  double bytes = 0.0;
  if (dflag) bytes += nmax * sizeof(double);
  if (eflag) bytes += nmax * sizeof(double);
  if (nvalues > 1) bytes += nmax*nvalues * sizeof(double);
  return bytes;
}

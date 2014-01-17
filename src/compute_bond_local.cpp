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

#define DELTA 10000

enum{DIST,ENG,FORCE};

/* ---------------------------------------------------------------------- */

ComputeBondLocal::ComputeBondLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute bond/local command");

  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Compute bond/local used when bonds are not allowed");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  bstyle = new int[nvalues];

  nvalues = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"dist") == 0) bstyle[nvalues++] = DIST;
    else if (strcmp(arg[iarg],"eng") == 0) bstyle[nvalues++] = ENG;
    else if (strcmp(arg[iarg],"force") == 0) bstyle[nvalues++] = FORCE;
    else error->all(FLERR,"Invalid keyword in compute bond/local command");
  }

  // set singleflag if need to call bond->single()

  singleflag = 0;
  for (int i = 0; i < nvalues; i++)
    if (bstyle[i] != DIST) singleflag = 1;

  nmax = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBondLocal::~ComputeBondLocal()
{
  memory->destroy(vector);
  memory->destroy(array);
  delete [] bstyle;
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::init()
{
  if (force->bond == NULL)
    error->all(FLERR,"No bond style is defined for compute bond/local");

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
}

/* ----------------------------------------------------------------------
   count bonds and compute bond info on this proc
   only count bond once if newton_bond is off
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if bond is deleted (type = 0), do not count
   if bond is turned off (type < 0), still count
   if flag is set, compute requested info about bond
   if bond is turned off (type < 0), energy = 0.0
------------------------------------------------------------------------- */

int ComputeBondLocal::compute_bonds(int flag)
{
  int i,m,n,atom1,atom2;
  double delx,dely,delz,rsq;
  double *dbuf,*ebuf,*fbuf;
  double *ptr;

  double **x = atom->x;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  Bond *bond = force->bond;
  double eng,fbond;

  m = n = 0;
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

        if (singleflag) {
          if (bond_type[atom1][i] > 0)
            eng = bond->single(bond_type[atom1][i],rsq,atom1,atom2,fbond);
          else eng = fbond = 0.0;
        }

        if (nvalues == 1) ptr = &vector[m];
        else ptr = array[m];

        for (n = 0; n < nvalues; n++) {
          switch (bstyle[n]) {
          case DIST:
            ptr[n] = sqrt(rsq);
            break;
          case ENG:
            ptr[n] = eng;
            break;
          case FORCE:
            ptr[n] = sqrt(rsq)*fbond;
            break;
          }
        }
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vector);
    memory->create(vector,nmax,"bond/local:vector");
    vector_local = vector;
  } else {
    memory->destroy(array);
    memory->create(array,nmax,nvalues,"bond/local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeBondLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}

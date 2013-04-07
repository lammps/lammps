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

#include "stdlib.h"
#include "string.h"
#include "fix_store.h"
#include "atom.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixStore::FixStore(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix store command");

  // syntax: id group style 0/1 nvalue

  restart_peratom = atoi(arg[3]);
  nvalues = atoi(arg[4]);

  vecflag = 0;
  if (nvalues == 1) vecflag = 1;

  // perform initial allocation of atom-based array
  // register with Atom class

  vstore = NULL;
  astore = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  if (restart_peratom) atom->add_callback(1);

  // zero the storage
  // since may be exchanged before filled by caller, e.g. fix store/force

  int nlocal = atom->nlocal;
  if (vecflag)
    for (int i = 0; i < nlocal; i++)
      vstore[i] = 0.0;
  else 
    for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < nvalues; j++)
	astore[i][j] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixStore::~FixStore()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  if (restart_peratom) atom->delete_callback(id,1);

  memory->destroy(vstore);
  memory->destroy(astore);
}

/* ---------------------------------------------------------------------- */

int FixStore::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixStore::memory_usage()
{
  double bytes = atom->nmax*nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixStore::grow_arrays(int nmax)
{
  if (vecflag) memory->grow(vstore,nmax,"store:vstore");
  else memory->grow(astore,nmax,nvalues,"store:astore");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixStore::copy_arrays(int i, int j, int delflag)
{
  if (vecflag) vstore[j] = vstore[i];
  else
    for (int m = 0; m < nvalues; m++)
      astore[j][m] = astore[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixStore::pack_exchange(int i, double *buf)
{
  if (vecflag) buf[0] = vstore[i];
  else
    for (int m = 0; m < nvalues; m++)
      buf[m] = astore[i][m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixStore::unpack_exchange(int nlocal, double *buf)
{
  if (vecflag) vstore[nlocal] = buf[0];
  else 
    for (int m = 0; m < nvalues; m++)
      astore[nlocal][m] = buf[m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixStore::pack_restart(int i, double *buf)
{
  buf[0] = nvalues+1;
  if (vecflag) buf[1] = vstore[i];
  else
    for (int m = 0; m < nvalues; m++)
      buf[m+1] = astore[i][m];
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixStore::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  if (vecflag) vstore[nlocal] = extra[nlocal][m];
  else
    for (int i = 0; i < nvalues; i++) 
      astore[nlocal][i] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixStore::maxsize_restart()
{
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixStore::size_restart(int nlocal)
{
  return nvalues+1;
}

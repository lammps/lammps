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

#include <stdlib.h>
#include <string.h>
#include "fix_store.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{GLOBAL,PERATOM};

/* ---------------------------------------------------------------------- */

FixStore::FixStore(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix store command");

  // 4th arg determines GLOBAL vs PERATOM values
  // syntax: id group style peratom 0/1 nvalue
  //   0/1 flag = not-store or store peratom values in restart file
  //   nvalue = # of peratom values, N=1 is vector, N>1 is array
  // syntax: id group style global nrow ncol
  //   Nrow by Ncol array of global values
  //   Ncol=1 is vector, Nrow>1 is array

  if (strcmp(arg[3],"peratom") == 0) flavor = PERATOM;
  else if (strcmp(arg[3],"global") == 0) flavor = GLOBAL;
  else error->all(FLERR,"Invalid fix store command");

  // GLOBAL values are always written to restart file
  // PERATOM restart_peratom is set by caller

  if (flavor == GLOBAL) {
    restart_global = 1;
    nrow = force->inumeric(FLERR,arg[4]);
    ncol = force->inumeric(FLERR,arg[5]);
    if (nrow <= 0 || ncol <= 0)
      error->all(FLERR,"Invalid fix store command");
    vecflag = 0;
    if (ncol == 1) vecflag = 1;
  } else {
    restart_peratom = force->inumeric(FLERR,arg[4]);
    nvalues = force->inumeric(FLERR,arg[5]);
    if (restart_peratom < 0 or restart_peratom > 1 || nvalues <= 0)
      error->all(FLERR,"Invalid fix store command");
    vecflag = 0;
    if (nvalues == 1) vecflag = 1;
  }

  vstore = NULL;
  astore = NULL;

  // allocate vector or array
  // for PERATOM, register with Atom class

  if (flavor == GLOBAL) {
    if (vecflag) memory->create(vstore,nrow,"fix/store:vstore");
    else memory->create(astore,nrow,ncol,"fix/store:astore");
  }
  if (flavor == PERATOM) {
    grow_arrays(atom->nmax);
    atom->add_callback(0);
    if (restart_peratom) atom->add_callback(1);
  }

  // zero the storage
  // PERATOM may be exchanged before filled by caller

  if (flavor == GLOBAL) {
    if (vecflag)
      for (int i = 0; i < nrow; i++) vstore[i] = 0.0;
    else
      for (int i = 0; i < nrow; i++)
        for (int j = 0; j < ncol; j++)
          astore[i][j] = 0.0;
  }
  if (flavor == PERATOM) {
    int nlocal = atom->nlocal;
    if (vecflag)
      for (int i = 0; i < nlocal; i++) vstore[i] = 0.0;
    else
      for (int i = 0; i < nlocal; i++)
        for (int j = 0; j < nvalues; j++)
          astore[i][j] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixStore::~FixStore()
{
  // unregister callbacks to this fix from Atom class

  if (flavor == PERATOM) {
    atom->delete_callback(id,0);
    if (restart_peratom) atom->delete_callback(id,1);
  }

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
   write global array to restart file
------------------------------------------------------------------------- */

void FixStore::write_restart(FILE *fp)
{
  int n = nrow*ncol;
  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    if (vecflag) fwrite(vstore,sizeof(double),n,fp);
    else fwrite(&astore[0][0],sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use global array from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixStore::restart(char *buf)
{
  // HOWTO insure size of buf is the same

  int n = nrow*ncol;
  double *dbuf = (double *) buf;
  if (vecflag) memcpy(vstore,dbuf,n*sizeof(double));
  else memcpy(&astore[0][0],dbuf,n*sizeof(double));
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

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixStore::memory_usage()
{
  double bytes;
  if (flavor == GLOBAL) bytes += nrow*ncol * sizeof(double);
  if (flavor == PERATOM) bytes += atom->nmax*nvalues * sizeof(double);
  return bytes;
}

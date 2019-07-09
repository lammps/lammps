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

#include <cstdlib>
#include <cstring>
#include "fix_store.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{UNKNOWN,GLOBAL,PERATOM};

/* ---------------------------------------------------------------------- */

FixStore::FixStore(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
vstore(NULL), astore(NULL), rbuf(NULL)
{
  if (narg != 6) error->all(FLERR,"Illegal fix store command");

  // 4th arg determines GLOBAL vs PERATOM values
  // syntax: id group style global nrow ncol
  //   Nrow by Ncol array of global values
  //   Ncol = 1 is vector, Ncol > 1 is array
  // syntax: id group style peratom 0/1 nvalues
  //   0/1 flag = not-store or store peratom values in restart file
  //   nvalues = # of peratom values, N = 1 is vector, N > 1 is array

  disable = 0;
  nvalues = vecflag = 0;
  flavor = UNKNOWN;

  if (strcmp(arg[3],"global") == 0) flavor = GLOBAL;
  else if (strcmp(arg[3],"peratom") == 0) flavor = PERATOM;
  else error->all(FLERR,"Illegal fix store command");

  // GLOBAL values are always written to restart file
  // PERATOM restart_peratom is set by caller

  if (flavor == GLOBAL) {
    restart_global = 1;
    nrow = force->inumeric(FLERR,arg[4]);
    ncol = force->inumeric(FLERR,arg[5]);
    if (nrow <= 0 || ncol <= 0)
      error->all(FLERR,"Illegal fix store command");
    vecflag = 0;
    if (ncol == 1) vecflag = 1;
  }
  if (flavor == PERATOM) {
    restart_peratom = force->inumeric(FLERR,arg[4]);
    nvalues = force->inumeric(FLERR,arg[5]);
    if (restart_peratom < 0 or restart_peratom > 1 || nvalues <= 0)
      error->all(FLERR,"Illegal fix store command");
    vecflag = 0;
    if (nvalues == 1) vecflag = 1;
  }

  vstore = NULL;
  astore = NULL;

  // allocate vector or array and restart buffer rbuf
  // for PERATOM, register with Atom class

  if (flavor == GLOBAL) {
    if (vecflag) memory->create(vstore,nrow,"fix/store:vstore");
    else memory->create(astore,nrow,ncol,"fix/store:astore");
    memory->create(rbuf,nrow*ncol+2,"fix/store:rbuf");
  }
  if (flavor == PERATOM) {
    grow_arrays(atom->nmax);
    atom->add_callback(0);
    if (restart_peratom) atom->add_callback(1);
    rbuf = NULL;
  }

  // zero the storage
  // PERATOM may be comm->exchanged before filled by caller

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
    maxexchange = nvalues;
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
  memory->destroy(rbuf);
}

/* ---------------------------------------------------------------------- */

int FixStore::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   reset size of global vector/array
   invoked by caller if size is unknown at time this fix is instantiated
   caller will do subsequent initialization
------------------------------------------------------------------------- */

void FixStore::reset_global(int nrow_caller, int ncol_caller)
{
  memory->destroy(vstore);
  memory->destroy(astore);
  memory->destroy(rbuf);
  vstore = NULL;
  astore = NULL;

  vecflag = 0;
  if (ncol_caller == 1) vecflag = 1;
  nrow = nrow_caller;
  ncol = ncol_caller;
  if (vecflag) memory->create(vstore,nrow,"fix/store:vstore");
  else memory->create(astore,nrow,ncol,"fix/store:astore");
  memory->create(rbuf,nrow*ncol+2,"fix/store:rbuf");
}

/* ----------------------------------------------------------------------
   write global array to restart file
------------------------------------------------------------------------- */

void FixStore::write_restart(FILE *fp)
{
  // fill rbuf with size and vec/array values

  rbuf[0] = nrow;
  rbuf[1] = ncol;
  if (vecflag) memcpy(&rbuf[2],vstore,nrow*sizeof(double));
  else memcpy(&rbuf[2],&astore[0][0],nrow*ncol*sizeof(double));

  int n = nrow*ncol + 2;
  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rbuf,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use global array from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixStore::restart(char *buf)
{
  // first 2 values in buf are vec/array sizes

  double *dbuf = (double *) buf;
  int nrow_restart = dbuf[0];
  int ncol_restart = dbuf[1];

  // if size of vec/array has changed,
  //   means the restart file is setting size of vec or array and doing init
  //   because caller did not know size at time this fix was instantiated
  // reallocate vstore or astore accordingly

  if (nrow != nrow_restart || ncol != ncol_restart) {
    memory->destroy(vstore);
    memory->destroy(astore);
    memory->destroy(rbuf);
    vstore = NULL;
    astore = NULL;

    vecflag = 0;
    if (ncol_restart == 1) vecflag = 1;
    nrow = nrow_restart;
    ncol = ncol_restart;
    if (vecflag) memory->create(vstore,nrow,"fix/store:vstore");
    else memory->create(astore,nrow,ncol,"fix/store:astore");
    memory->create(rbuf,nrow*ncol+2,"fix/store:rbuf");
  }

  int n = nrow*ncol;
  if (vecflag) memcpy(vstore,&dbuf[2],n*sizeof(double));
  else memcpy(&astore[0][0],&dbuf[2],n*sizeof(double));
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

void FixStore::copy_arrays(int i, int j, int /*delflag*/)
{
  if (disable) return;

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
  if (disable) return 0;

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
  if (disable) return 0;

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
  if (disable) {
    buf[0] = 0;
    return 1;
  }

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
  if (disable) return;

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
  if (disable) return 1;
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixStore::size_restart(int /*nlocal*/)
{
  if (disable) return 1;
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixStore::memory_usage()
{
  double bytes = 0.0;
  if (flavor == GLOBAL) bytes += nrow*ncol * sizeof(double);
  if (flavor == PERATOM) bytes += atom->nmax*nvalues * sizeof(double);
  return bytes;
}

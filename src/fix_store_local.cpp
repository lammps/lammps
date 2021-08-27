/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_store_local.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 1000

/* ---------------------------------------------------------------------- */

FixStoreLocal::FixStoreLocal(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nvalues(0), vector(nullptr), array(nullptr)
{
  if (narg != 5) error->all(FLERR, "Illegal fix store/local command");
  local_flag = 1;

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, "Illegal fix store/local command");
  local_freq = nevery;

  nvalues = utils::inumeric(FLERR, arg[4], false, lmp);
  
  if (nvalues <= 0) error->all(FLERR, "Illegal fix store/local command");
  if (nvalues == 1)
    size_local_cols = 0;
  else
    size_local_cols = nvalues;
  size_local_rows = 0;
  
  vector = nullptr;
  array = nullptr;  
  nmax = 0;
  ncount = 0;
}

/* ---------------------------------------------------------------------- */

FixStoreLocal::~FixStoreLocal()
{
  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixStoreLocal::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStoreLocal::add_data(double *input_data, int i, int j)
{
  int *mask = atom->mask;
  if (!(mask[i] & groupbit)) return;
  if (!(mask[j] & groupbit)) return;
  
  if (ncount >= nmax) reallocate(ncount);
  
  // fill vector or array with local values
  if (nvalues == 1) {
    vector[ncount] = input_data[0];
  } else {
    for (int i = 0; i < nvalues; i++) array[ncount][i] = input_data[i];
  }

  ncount += 1;
}

/* ---------------------------------------------------------------------- */

void FixStoreLocal::post_force(int /*vflag*/)
{
  if (update->ntimestep % nevery == 0) {
    size_local_rows = ncount;
    ncount = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreLocal::reallocate(int n)
{
  // grow vector or array
  while (nmax <= n) nmax += DELTA;

  if (nvalues == 1) {
    memory->grow(vector, nmax, "fix_store_local:vector");
    vector_local = vector;
  } else {
    memory->grow(array, nmax, nvalues, "fix_store_local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   write global array to restart file 
------------------------------------------------------------------------- */

void FixStoreLocal::write_restart(FILE *fp)
{
  // fill rbuf with size and vec/array values

  rbuf[0] = nmax;
  rbuf[1] = nvalues;
  if (nvalues == 1) memcpy(&rbuf[2],vector,sizeof(double)*nmax);
  else memcpy(&rbuf[2],&array[0][0],sizeof(double)*nmax*nvalues);

  int n = nmax*nvalues + 2;
  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rbuf,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use global array from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixStoreLocal::restart(char *buf)
{
  // first 2 values in buf are vec/array sizes

  double *dbuf = (double *) buf;
  int nrow_restart = dbuf[0];
  int ncol_restart = dbuf[1];

  // if size of vec/array has changed,
  //   means the restart file is setting size of vec or array and doing init
  //   because caller did not know size at time this fix was instantiated
  // reallocate vector or array accordingly

  if (nmax != nrow_restart || nvalues != ncol_restart) {
    memory->destroy(vector);
    memory->destroy(array);
    memory->destroy(rbuf);
    vector = nullptr;
    array = nullptr;

    nmax = nrow_restart;
    nvalues = ncol_restart;
    if (nvalues == 1) memory->create(vector,nmax,"fix/store/local:vector");
    else memory->create(array,nmax,nvalues,"fix/store/local:array");
    memory->create(rbuf,nmax*nvalues+2,"fix/store:rbuf");
  }

  int n = nmax*nvalues;
  if (nvalues == 1) memcpy(vector,&dbuf[2],n*sizeof(double));
  else memcpy(&array[0][0],&dbuf[2],n*sizeof(double));
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixStoreLocal::maxsize_restart()
{
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixStoreLocal::size_restart(int /*nlocal*/)
{
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double FixStoreLocal::memory_usage()
{
  double bytes = (double) nmax * (double) nvalues * sizeof(double);
  return bytes;
}


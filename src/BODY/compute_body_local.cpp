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

#include "compute_body_local.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "atom_vec_body.h"
#include "body.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

enum{ID,TYPE,INDEX};

/* ---------------------------------------------------------------------- */

ComputeBodyLocal::ComputeBodyLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), which(NULL), index(NULL), avec(NULL), bptr(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute body/local command");

  local_flag = 1;
  nvalues = narg - 3;

  which = new int[nvalues];
  index = new int[nvalues];
  nvalues = 0;

  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"id") == 0) which[nvalues++] = ID;
    else if (strcmp(arg[iarg],"type") == 0) which[nvalues++] = TYPE;
    else {
      which[nvalues] = INDEX;
      index[nvalues] = force->inumeric(FLERR,arg[iarg]) - 1;
      nvalues++;
    }
  }

  avec = (AtomVecBody *) atom->style_match("body");
  if (!avec) error->all(FLERR,"Compute body/local requires atom style body");
  bptr = avec->bptr;

  int indexmax = bptr->noutcol();
  for (int i = 0; i < nvalues; i++) {
    if (which[i] == INDEX && (index[i] < 0 || index[i] >= indexmax))
      error->all(FLERR,"Invalid index in compute body/local command");
  }

  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  nmax = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBodyLocal::~ComputeBodyLocal()
{
  delete [] which;
  delete [] index;
  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeBodyLocal::init()
{
  // if non-body particles in group insure only indices 1,2,3 are used

  int nonbody = 0;
  int *mask = atom->mask;
  int *body = atom->body;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (body[i] < 0) nonbody = 1;

  int flag;
  MPI_Allreduce(&nonbody,&flag,1,MPI_INT,MPI_SUM,world);

  if (flag) {
    for (int i = 0; i < nvalues; i++)
      if (which[i] == INDEX && index[i] > 2)
        error->all(FLERR,"Invalid index for non-body particles "
                   "in compute body/local command");
  }

  // do initial memory allocation so that memory_usage() is correct

  int ncount = compute_body(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeBodyLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute body info

  int ncount = compute_body(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_body(1);
}

/* ----------------------------------------------------------------------
   count body info and compute body info on this proc
   flag = 0, only count
   flag = 1, also compute
------------------------------------------------------------------------- */

int ComputeBodyLocal::compute_body(int flag)
{
  // perform count

  int *mask = atom->mask;
  int *body = atom->body;
  int nlocal = atom->nlocal;

  int ncount = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (body[i] < 0) ncount++;
      else ncount += bptr->noutrow(body[i]);
    }

  if (flag == 0) return ncount;

  // perform computation and fill output vector/array

  int m,n,ibonus;
  double *values = new double[bptr->noutcol()];

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *type = atom->type;

  ncount = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (body[i] < 0) {
        if (nvalues == 1) {
          if (which[0] == ID) vector[ncount] = tag[i];
          else if (which[0] == TYPE) vector[ncount] = type[i];
          else vector[ncount] = x[i][index[0]];
        } else {
          for (m = 0; m < nvalues; m++) {
            if (which[m] == ID) array[ncount][m] = tag[i];
            else if (which[m] == TYPE) array[ncount][m] = type[i];
            else array[ncount][m] = x[i][index[m]];
          }
        }
        ncount++;

      } else {
        ibonus = body[i];
        n = bptr->noutrow(ibonus);
        for (int j = 0; j < n; j++) {
          bptr->output(ibonus,j,values);
          if (nvalues == 1) {
            if (which[0] == ID) vector[ncount] = tag[i];
            else if (which[0] == TYPE) vector[ncount] = type[i];
            else vector[ncount] = values[index[0]];
          } else {
            for (m = 0; m < nvalues; m++) {
              if (which[m] == ID) array[ncount][m] = tag[i];
              else if (which[m] == TYPE) array[ncount][m] = type[i];
              else array[ncount][m] = values[index[m]];
            }
          }
          ncount++;
        }
      }
    }
  }

  delete [] values;
  return ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeBodyLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vector);
    memory->create(vector,nmax,"body/local:vector");
    vector_local = vector;
  } else {
    memory->destroy(array);
    memory->create(array,nmax,nvalues,"body/local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeBodyLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}

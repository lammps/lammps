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

#include <cstring>
#include <cstdlib>
#include "compute_chunk_spread_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "compute_chunk_atom.h"
#include "input.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{COMPUTE,FIX};

#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

ComputeChunkSpreadAtom::
ComputeChunkSpreadAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idchunk(NULL), ids(NULL), which(NULL), argindex(NULL), value2index(NULL)
{
  if (narg < 5) error->all(FLERR,"Illegal compute chunk/spread/atom command");

  // ID of compute chunk/atom

  int n = strlen(arg[3]) + 1;
  idchunk = new char[n];
  strcpy(idchunk,arg[3]);
  init_chunk();

  // expand args if any have wildcard character "*"

  int iarg = 4;
  int expand = 0;
  char **earg;
  int nargnew = input->expand_args(narg-iarg,&arg[iarg],1,earg);

  if (earg != &arg[iarg]) expand = 1;
  arg = earg;

  // parse values

  which = new int[nargnew];
  argindex = new int[nargnew];
  ids = new char*[nargnew];
  value2index = new int[nargnew];
  nvalues = 0;

  iarg = 0;
  while (iarg < nargnew) {
    ids[nvalues] = NULL;

    if (strncmp(arg[iarg],"c_",2) == 0 ||
        strncmp(arg[iarg],"f_",2) == 0) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal compute chunk/spread/atom command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else error->all(FLERR,"Illegal compute chunk/spread/atom command");

    iarg++;
  }

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // setup and error check
  // for compute, must calculate per-chunk values, i.e. style ends in "/chunk"
  // for fix, assume a global vector or array is per-chunk data

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute chunk/spread/atom "
                   "does not exist");

      char *ptr = strstr(modify->compute[icompute]->style,"/chunk");
      if (!ptr || (ptr !=  modify->compute[icompute]->style +
                   strlen(modify->compute[icompute]->style) - strlen("/chunk")))
        error->all(FLERR,"Compute for compute chunk/spread/atom "
                   "does not calculate per-chunk values");

      if (argindex[i] == 0) {
        if (!modify->compute[icompute]->vector_flag)
          error->all(FLERR,"Compute chunk/spread/atom compute "
                     "does not calculate global vector");
      } else {
        if (!modify->compute[icompute]->array_flag)
          error->all(FLERR,"Compute chunk/spread/atom compute "
                     "does not calculate global array");
        if (argindex[i] > modify->compute[icompute]->size_array_cols)
          error->all(FLERR,"Compute chunk/spread/atom compute array "
                     "is accessed out-of-range");
      }

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute chunk/spread/atom does not exist");
      if (argindex[i] == 0) {
        if (!modify->fix[ifix]->vector_flag)
          error->all(FLERR,"Compute chunk/spread/atom fix "
                     "does not calculate global vector");
      } else {
        if (!modify->fix[ifix]->array_flag)
          error->all(FLERR,"Compute chunk/spread/atom fix "
                     "does not calculate global array");
        if (argindex[i] > modify->fix[ifix]->size_array_cols)
          error->all(FLERR,"Compute chunk/spread/atom fix array "
                     "is accessed out-of-range");
      }
    }
  }

  // this compute produces a peratom vector or array

  peratom_flag = 1;
  if (nvalues == 1) size_peratom_cols = 0;
  else size_peratom_cols = nvalues;

  // per-atom vector or array

  nmax = 0;
  vector_atom = NULL;
  array_atom = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeChunkSpreadAtom::~ComputeChunkSpreadAtom()
{
  delete [] idchunk;

  delete [] which;
  delete [] argindex;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;
  delete [] value2index;

  memory->destroy(vector_atom);
  memory->destroy(array_atom);
}

/* ---------------------------------------------------------------------- */

void ComputeChunkSpreadAtom::init()
{
  init_chunk();

  // set indices of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute chunk/spread/atom "
                   "does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute chunk/spread/atom does not exist");
      value2index[m] = ifix;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeChunkSpreadAtom::init_chunk()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute < 0)
    error->all(FLERR,"Chunk/atom compute does not exist for compute chunk/spread/atom");
  cchunk = (ComputeChunkAtom *) modify->compute[icompute];
  if (strcmp(cchunk->style,"chunk/atom") != 0)
    error->all(FLERR,"Compute chunk/spread/atom does not use chunk/atom compute");
}

/* ---------------------------------------------------------------------- */

void ComputeChunkSpreadAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow local vector_atom or array_atom if necessary

  if (atom->nmax > nmax) {
    if (nvalues == 1) {
      memory->destroy(vector_atom);
      nmax = atom->nmax;
      memory->create(vector_atom,nmax,"chunk/spread/atom:vector_atom");
    } else {
      memory->destroy(array_atom);
      nmax = atom->nmax;
      memory->create(array_atom,nmax,nvalues,"chunk/spread/atom:array_atom");
    }
  }

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  int nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;

  // loop over values, access compute or fix
  // loop over atoms, use chunk ID of each atom to store value from compute/fix

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int i,m,n,index,nstride;
  double *ptr;

  for (m = 0; m < nvalues; m++) {
    n = value2index[m];

    // copy compute/fix values into vector_atom or array_atom
    // nstride between values for each atom

    if (nvalues == 1) {
      ptr = vector_atom;
      nstride = 1;
    } else {
      ptr = &array_atom[0][m];
      nstride = nvalues;
    }

    // invoke compute if not previously invoked

    if (which[m] == COMPUTE) {
      Compute *compute = modify->compute[n];

      if (argindex[m] == 0) {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        double *cvector = compute->vector;
        for (i = 0; i < nlocal; i++, ptr += nstride) {
          *ptr = 0.0;
          if (!(mask[i] & groupbit)) continue;
          index = ichunk[i]-1;
          if (index < 0 || index >= nchunk) continue;
          *ptr = cvector[index];
        }

      } else {
        if (!(compute->invoked_flag & INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= INVOKED_ARRAY;
        }
        int icol = argindex[m]-1;
        double **carray = compute->array;
        for (i = 0; i < nlocal; i++, ptr += nstride) {
          *ptr = 0.0;
          if (!(mask[i] & groupbit)) continue;
          index = ichunk[i]-1;
          if (index < 0 || index >= nchunk) continue;
          *ptr = carray[index][icol];
        }
      }

    // access fix data, check if fix frequency is a match
    // are assuming the fix global vector/array is per-chunk data
    // check if index exceeds fix output length/rows

    } else if (which[m] == FIX) {
      Fix *fix = modify->fix[n];
      if (update->ntimestep % fix->global_freq)
        error->all(FLERR,"Fix used in compute chunk/spread/atom not "
                   "computed at compatible time");

      if (argindex[m] == 0) {
        int nfix = fix->size_vector;
        for (i = 0; i < nlocal; i++, ptr += nstride) {
          *ptr = 0.0;
          if (!(mask[i] & groupbit)) continue;
          index = ichunk[i]-1;
          if (index < 0 || index >= nchunk || index >= nfix) continue;
          *ptr = fix->compute_vector(index);
        }

      } else {
        int icol = argindex[m]-1;
        int nfix = fix->size_array_rows;
        for (i = 0; i < nlocal; i++, ptr += nstride) {
          *ptr = 0.0;
          if (!(mask[i] & groupbit)) continue;
          index = ichunk[i]-1;
          if (index < 0 || index >= nchunk || index >= nfix) continue;
          *ptr = fix->compute_array(index,icol);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeChunkSpreadAtom::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}

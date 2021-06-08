// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_chunk_spread_atom.h"

#include "arg_info.h"
#include "atom.h"
#include "compute.h"
#include "compute_chunk_atom.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeChunkSpreadAtom::
ComputeChunkSpreadAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idchunk(nullptr), ids(nullptr), which(nullptr), argindex(nullptr), value2index(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal compute chunk/spread/atom command");

  // ID of compute chunk/atom

  idchunk = utils::strdup(arg[3]);
  init_chunk();

  // expand args if any have wildcard character "*"

  int iarg = 4;
  int expand = 0;
  char **earg;
  int nargnew = utils::expand_args(FLERR,narg-iarg,&arg[iarg],1,earg,lmp);

  if (earg != &arg[iarg]) expand = 1;
  arg = earg;

  // parse values

  which = new int[nargnew];
  argindex = new int[nargnew];
  ids = new char*[nargnew];
  value2index = new int[nargnew];
  nvalues = 0;

  for (iarg = 0; iarg < nargnew; iarg++) {
    ids[nvalues] = nullptr;

    ArgInfo argi(arg[iarg], ArgInfo::COMPUTE|ArgInfo::FIX);

    which[nvalues] = argi.get_type();
    argindex[nvalues] = argi.get_index1();
    ids[nvalues] = argi.copy_name();

    if ((which[nvalues] == ArgInfo::UNKNOWN) || (which[nvalues] == ArgInfo::NONE)
        || (argi.get_dim() > 1))
      error->all(FLERR,"Illegal compute chunk/spread/atom command");

    nvalues++;
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
    if (which[i] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute chunk/spread/atom "
                   "does not exist");

      if (!utils::strmatch(modify->compute[icompute]->style,"/chunk$"))
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

    } else if (which[i] == ArgInfo::FIX) {
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
  vector_atom = nullptr;
  array_atom = nullptr;
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
    if (which[m] == ArgInfo::COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute chunk/spread/atom "
                   "does not exist");
      value2index[m] = icompute;

    } else if (which[m] == ArgInfo::FIX) {
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

    if (which[m] == ArgInfo::COMPUTE) {
      Compute *compute = modify->compute[n];

      if (argindex[m] == 0) {
        if (!(compute->invoked_flag & Compute::INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= Compute::INVOKED_VECTOR;
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
        if (!(compute->invoked_flag & Compute::INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= Compute::INVOKED_ARRAY;
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

    } else if (which[m] == ArgInfo::FIX) {
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
  double bytes = (double)nmax*nvalues * sizeof(double);
  return bytes;
}

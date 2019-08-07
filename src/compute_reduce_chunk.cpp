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

#include "compute_reduce_chunk.h"
#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "compute_chunk_atom.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{SUM,MINN,MAXX};
enum{UNKNOWN=-1,COMPUTE,FIX,VARIABLE};

#define INVOKED_PERATOM 8

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeReduceChunk::ComputeReduceChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  which(NULL), argindex(NULL), value2index(NULL), idchunk(NULL), ids(NULL),
  vlocal(NULL), vglobal(NULL), alocal(NULL), aglobal(NULL), varatom(NULL)
{
  if (narg < 6) error->all(FLERR,"Illegal compute reduce/chunk command");

  // ID of compute chunk/atom

  int n = strlen(arg[3]) + 1;
  idchunk = new char[n];
  strcpy(idchunk,arg[3]);
  init_chunk();

  // mode

  if (strcmp(arg[4],"sum") == 0) mode = SUM;
  else if (strcmp(arg[4],"min") == 0) mode = MINN;
  else if (strcmp(arg[4],"max") == 0) mode = MAXX;
  else error->all(FLERR,"Illegal compute reduce/chunk command");

  int iarg = 5;

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  int nargnew = input->expand_args(narg-iarg,&arg[iarg],1,earg);

  if (earg != &arg[iarg]) expand = 1;
  arg = earg;

  // parse values until one isn't recognized

  which = new int[nargnew];
  argindex = new int[nargnew];
  ids = new char*[nargnew];
  value2index = new int[nargnew];
  for (int i=0; i < nargnew; ++i) {
    which[i] = argindex[i] = value2index[i] = UNKNOWN;
    ids[i] = NULL;
  }
  nvalues = 0;

  iarg = 0;
  while (iarg < nargnew) {
    ids[nvalues] = NULL;

    if (strncmp(arg[iarg],"c_",2) == 0 ||
               strncmp(arg[iarg],"f_",2) == 0 ||
               strncmp(arg[iarg],"v_",2) == 0) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal compute reduce/chunk command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else error->all(FLERR,"Illegal compute reduce/chunk command");

    iarg++;
  }

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // error check

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute reduce/chunk does not exist");
      if (!modify->compute[icompute]->peratom_flag)
        error->all(FLERR,"Compute reduce/chunk compute does not "
                   "calculate per-atom values");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_peratom_cols != 0)
        error->all(FLERR,"Compute reduce/chunk compute does not "
                   "calculate a per-atom vector");
      if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
        error->all(FLERR,"Compute reduce/chunk compute does not "
                   "calculate a per-atom array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_peratom_cols)
        error->all(FLERR,
                   "Compute reduce/chunk compute array is accessed out-of-range");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute reduce/chunk does not exist");
      if (!modify->fix[ifix]->peratom_flag)
        error->all(FLERR,"Compute reduce/chunk fix does not "
                   "calculate per-atom values");
      if (argindex[i] == 0 &&
          modify->fix[ifix]->size_peratom_cols != 0)
        error->all(FLERR,"Compute reduce/chunk fix does not "
                   "calculate a per-atom vector");
      if (argindex[i] && modify->fix[ifix]->size_peratom_cols == 0)
        error->all(FLERR,"Compute reduce/chunk fix does not "
                   "calculate a per-atom array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_peratom_cols)
        error->all(FLERR,"Compute reduce/chunk fix array is "
                   "accessed out-of-range");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for compute reduce/chunk does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
        error->all(FLERR,"Compute reduce/chunk variable is "
                   "not atom-style variable");
    }
  }

  // this compute produces either a vector or array

  if (nvalues == 1) {
    vector_flag = 1;
    size_vector_variable = 1;
    extvector = 0;
  } else {
    array_flag = 1;
    size_array_rows_variable = 1;
    size_array_cols = nvalues;
    extarray = 0;
  }

  // setup

  if (mode == SUM) initvalue = 0.0;
  else if (mode == MINN) initvalue = BIG;
  else if (mode == MAXX) initvalue = -BIG;

  maxchunk = 0;
  vlocal = vglobal = NULL;
  alocal = aglobal = NULL;

  maxatom = 0;
  varatom = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeReduceChunk::~ComputeReduceChunk()
{
  delete [] idchunk;

  delete [] which;
  delete [] argindex;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;

  memory->destroy(vlocal);
  memory->destroy(vglobal);
  memory->destroy(alocal);
  memory->destroy(aglobal);

  memory->destroy(varatom);
}

/* ---------------------------------------------------------------------- */

void ComputeReduceChunk::init()
{
  init_chunk();

  // set indices of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute reduce/chunk does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute reduce/chunk does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for compute reduce/chunk does not exist");
      value2index[m] = ivariable;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeReduceChunk::init_chunk()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute < 0)
    error->all(FLERR,"Chunk/atom compute does not exist for "
               "compute reduce/chunk");
  cchunk = (ComputeChunkAtom *) modify->compute[icompute];
  if (strcmp(cchunk->style,"chunk/atom") != 0)
    error->all(FLERR,"Compute reduce/chunk does not use chunk/atom compute");
}

/* ---------------------------------------------------------------------- */

void ComputeReduceChunk::compute_vector()
{
  invoked_vector = update->ntimestep;

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  ichunk = cchunk->ichunk;
  if (!nchunk) return;

  size_vector = nchunk;

  if (nchunk > maxchunk) {
    memory->destroy(vlocal);
    memory->destroy(vglobal);
    maxchunk = nchunk;
    memory->create(vlocal,maxchunk,"reduce/chunk:vlocal");
    memory->create(vglobal,maxchunk,"reduce/chunk:vglobal");
    vector = vglobal;
  }

  // perform local reduction of single peratom value

  compute_one(0,vlocal,1);

  // reduce the per-chunk values across all procs

  if (mode == SUM)
    MPI_Allreduce(vlocal,vglobal,nchunk,MPI_DOUBLE,MPI_SUM,world);
  else if (mode == MINN)
    MPI_Allreduce(vlocal,vglobal,nchunk,MPI_DOUBLE,MPI_MIN,world);
  else if (mode == MAXX)
    MPI_Allreduce(vlocal,vglobal,nchunk,MPI_DOUBLE,MPI_MAX,world);
}

/* ---------------------------------------------------------------------- */

void ComputeReduceChunk::compute_array()
{
  invoked_array = update->ntimestep;

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  ichunk = cchunk->ichunk;
  if (!nchunk) return;

  size_array_rows = nchunk;

  if (nchunk > maxchunk) {
    memory->destroy(alocal);
    memory->destroy(aglobal);
    maxchunk = nchunk;
    memory->create(alocal,maxchunk,nvalues,"reduce/chunk:alocal");
    memory->create(aglobal,maxchunk,nvalues,"reduce/chunk:aglobal");
    array = aglobal;
  }

  // perform local reduction of all peratom values

  for (int m = 0; m < nvalues; m++) compute_one(m,&alocal[0][m],nvalues);

  // reduce the per-chunk values across all procs

  if (mode == SUM)
    MPI_Allreduce(&alocal[0][0],&aglobal[0][0],nchunk*nvalues,
                  MPI_DOUBLE,MPI_SUM,world);
  else if (mode == MINN)
    MPI_Allreduce(&alocal[0][0],&aglobal[0][0],nchunk*nvalues,
                  MPI_DOUBLE,MPI_MIN,world);
  else if (mode == MAXX)
    MPI_Allreduce(&alocal[0][0],&aglobal[0][0],nchunk*nvalues,
                  MPI_DOUBLE,MPI_MAX,world);
}

/* ---------------------------------------------------------------------- */

void ComputeReduceChunk::compute_one(int m, double *vchunk, int nstride)
{
  // initialize per-chunk values in accumulation vector

  for (int i = 0; i < nchunk; i += nstride) vchunk[i] = initvalue;

  // loop over my atoms
  // use peratom input and chunk ID of each atom to update vector

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int index = -1;
  int vidx = value2index[m];

  // initialization in case it has not yet been run, e.g. when
  // the compute was invoked right after it has been created

  if (vidx == UNKNOWN) {
    init();
    vidx = value2index[m];
  }

  if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[vidx];

    if (!(compute->invoked_flag & INVOKED_PERATOM)) {
      compute->compute_peratom();
      compute->invoked_flag |= INVOKED_PERATOM;
    }

    if (argindex[m] == 0) {
      double *vcompute = compute->vector_atom;
      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        index = ichunk[i]-1;
        if (index < 0) continue;
        combine(vchunk[index*nstride],vcompute[i]);
      }
    } else {
      double **acompute = compute->array_atom;
      int argindexm1 = argindex[m] - 1;
      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        index = ichunk[i]-1;
        if (index < 0) continue;
        combine(vchunk[index*nstride],acompute[i][argindexm1]);
      }
    }

  // access fix fields, check if fix frequency is a match

  } else if (which[m] == FIX) {
    Fix *fix = modify->fix[vidx];
    if (update->ntimestep % fix->peratom_freq)
      error->all(FLERR,"Fix used in compute reduce/chunk not "
                 "computed at compatible time");

    if (argindex[m] == 0) {
      double *vfix = fix->vector_atom;
      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        index = ichunk[i]-1;
        if (index < 0) continue;
        combine(vchunk[index*nstride],vfix[i]);
      }
    } else {
      double **afix = fix->array_atom;
      int argindexm1 = argindex[m] - 1;
      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        index = ichunk[i]-1;
        if (index < 0) continue;
        combine(vchunk[index*nstride],afix[i][argindexm1]);
      }
    }

  // evaluate atom-style variable

  } else if (which[m] == VARIABLE) {
    if (atom->nmax > maxatom) {
      memory->destroy(varatom);
      maxatom = atom->nmax;
      memory->create(varatom,maxatom,"reduce/chunk:varatom");
    }

    input->variable->compute_atom(vidx,igroup,varatom,1,0);
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      index = ichunk[i]-1;
      if (index < 0) continue;
      combine(vchunk[index*nstride],varatom[i]);
    }
  }
}

/* ----------------------------------------------------------------------
   combine two values according to reduction mode
------------------------------------------------------------------------- */

void ComputeReduceChunk::combine(double &one, double two)
{
  if (mode == SUM) one += two;
  else if (mode == MINN) {
    if (two < one) one = two;
  } else if (mode == MAXX) {
    if (two > one) one = two;
  }
}

/* ----------------------------------------------------------------------
   lock methods: called by fix ave/time
   these methods insure vector/array size is locked for Nfreq epoch
     by passing lock info along to compute chunk/atom
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   increment lock counter
------------------------------------------------------------------------- */

void ComputeReduceChunk::lock_enable()
{
  cchunk->lockcount++;
}

/* ----------------------------------------------------------------------
   decrement lock counter in compute chunk/atom, if it still exists
------------------------------------------------------------------------- */

void ComputeReduceChunk::lock_disable()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute >= 0) {
    cchunk = (ComputeChunkAtom *) modify->compute[icompute];
    cchunk->lockcount--;
  }
}

/* ----------------------------------------------------------------------
   calculate and return # of chunks = length of vector/array
------------------------------------------------------------------------- */

int ComputeReduceChunk::lock_length()
{
  nchunk = cchunk->setup_chunks();
  return nchunk;
}

/* ----------------------------------------------------------------------
   set the lock from startstep to stopstep
------------------------------------------------------------------------- */

void ComputeReduceChunk::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  cchunk->lock(fixptr,startstep,stopstep);
}

/* ----------------------------------------------------------------------
   unset the lock
------------------------------------------------------------------------- */

void ComputeReduceChunk::unlock(Fix *fixptr)
{
  cchunk->unlock(fixptr);
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeReduceChunk::memory_usage()
{
  double bytes = (bigint) maxatom * sizeof(double);
  if (nvalues == 1) bytes += (bigint) maxchunk * 2 * sizeof(double);
  else bytes += (bigint) maxchunk * nvalues * 2 * sizeof(double);
  return bytes;
}

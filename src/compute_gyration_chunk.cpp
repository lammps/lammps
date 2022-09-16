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

#include "compute_gyration_chunk.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGyrationChunk::ComputeGyrationChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idchunk(nullptr), massproc(nullptr), masstotal(nullptr), com(nullptr), comall(nullptr),
  rg(nullptr), rgall(nullptr), rgt(nullptr), rgtall(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute gyration/chunk command");

  // ID of compute chunk/atom

  idchunk = utils::strdup(arg[3]);

  ComputeGyrationChunk::init();

  // optional args

  tensor = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"tensor") == 0) {
      tensor = 1;
      iarg++;
    } else error->all(FLERR,"Illegal compute gyration/chunk command");
  }

  if (tensor) {
    array_flag = 1;
    size_array_cols = 6;
    size_array_rows = 0;
    size_array_rows_variable = 1;
    extarray = 0;
  } else {
    vector_flag = 1;
    size_vector = 0;
    size_vector_variable = 1;
    extvector = 0;
  }

  // chunk-based data

  nchunk = 1;
  maxchunk = 0;
  allocate();
}

/* ---------------------------------------------------------------------- */

ComputeGyrationChunk::~ComputeGyrationChunk()
{
  delete [] idchunk;
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(rg);
  memory->destroy(rgall);
  memory->destroy(rgt);
  memory->destroy(rgtall);
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationChunk::init()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute < 0)
    error->all(FLERR,"Chunk/atom compute does not exist for "
               "compute gyration/chunk");
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->compute[icompute]);
  if (strcmp(cchunk->style,"chunk/atom") != 0)
    error->all(FLERR,"Compute gyration/chunk does not use chunk/atom compute");
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationChunk::compute_vector()
{
  int i,index;
  double dx,dy,dz,massone;
  double unwrap[3];

  invoked_array = update->ntimestep;

  com_chunk();
  int *ichunk = cchunk->ichunk;

  for (i = 0; i < nchunk; i++) rg[i] = 0.0;

  // compute Rg for each chunk

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i]-1;
      if (index < 0) continue;
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      rg[index] += (dx*dx + dy*dy + dz*dz) * massone;
    }

  MPI_Allreduce(rg,rgall,nchunk,MPI_DOUBLE,MPI_SUM,world);

  for (i = 0; i < nchunk; i++)
    if (masstotal[i] > 0.0)
      rgall[i] = sqrt(rgall[i]/masstotal[i]);
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationChunk::compute_array()
{
  int i,j,index;
  double dx,dy,dz,massone;
  double unwrap[3];

  invoked_array = update->ntimestep;

  com_chunk();
  int *ichunk = cchunk->ichunk;

  for (i = 0; i < nchunk; i++)
    for (j = 0; j < 6; j++) rgt[i][j] = 0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i]-1;
      if (index < 0) continue;
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      rgt[index][0] += dx*dx * massone;
      rgt[index][1] += dy*dy * massone;
      rgt[index][2] += dz*dz * massone;
      rgt[index][3] += dx*dy * massone;
      rgt[index][4] += dx*dz * massone;
      rgt[index][5] += dy*dz * massone;
    }

  if (nchunk)
    MPI_Allreduce(&rgt[0][0],&rgtall[0][0],nchunk*6,MPI_DOUBLE,MPI_SUM,world);

  for (i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      for (j = 0; j < 6; j++)
        rgtall[i][j] = rgtall[i][j]/masstotal[i];
    }
  }
}


/* ----------------------------------------------------------------------
   calculate per-chunk COM, used by both scalar and tensor
------------------------------------------------------------------------- */

void ComputeGyrationChunk::com_chunk()
{
  int index;
  double massone;
  double unwrap[3];

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;

  if (nchunk > maxchunk) allocate();
  if (tensor) size_array_rows = nchunk;
  else size_vector = nchunk;

  // zero local per-chunk values

  for (int i = 0; i < nchunk; i++) {
    massproc[i] = 0.0;
    com[i][0] = com[i][1] = com[i][2] = 0.0;
  }

  // compute COM for each chunk

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i]-1;
      if (index < 0) continue;
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      domain->unmap(x[i],image[i],unwrap);
      massproc[index] += massone;
      com[index][0] += unwrap[0] * massone;
      com[index][1] += unwrap[1] * massone;
      com[index][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(massproc,masstotal,nchunk,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&com[0][0],&comall[0][0],3*nchunk,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      comall[i][0] /= masstotal[i];
      comall[i][1] /= masstotal[i];
      comall[i][2] /= masstotal[i];
    }
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

void ComputeGyrationChunk::lock_enable()
{
  cchunk->lockcount++;
}

/* ----------------------------------------------------------------------
   decrement lock counter in compute chunk/atom, it if still exists
------------------------------------------------------------------------- */

void ComputeGyrationChunk::lock_disable()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute >= 0) {
    cchunk = dynamic_cast<ComputeChunkAtom *>(modify->compute[icompute]);
    cchunk->lockcount--;
  }
}

/* ----------------------------------------------------------------------
   calculate and return # of chunks = length of vector/array
------------------------------------------------------------------------- */

int ComputeGyrationChunk::lock_length()
{
  nchunk = cchunk->setup_chunks();
  return nchunk;
}

/* ----------------------------------------------------------------------
   set the lock from startstep to stopstep
------------------------------------------------------------------------- */

void ComputeGyrationChunk::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  cchunk->lock(fixptr,startstep,stopstep);
}

/* ----------------------------------------------------------------------
   unset the lock
------------------------------------------------------------------------- */

void ComputeGyrationChunk::unlock(Fix *fixptr)
{
  cchunk->unlock(fixptr);
}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeGyrationChunk::allocate()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(rg);
  memory->destroy(rgall);
  memory->destroy(rgt);
  memory->destroy(rgtall);
  maxchunk = nchunk;
  memory->create(massproc,maxchunk,"gyration/chunk:massproc");
  memory->create(masstotal,maxchunk,"gyration/chunk:masstotal");
  memory->create(com,maxchunk,3,"gyration/chunk:com");
  memory->create(comall,maxchunk,3,"gyration/chunk:comall");
  if (tensor) {
    memory->create(rgt,maxchunk,6,"gyration/chunk:rgt");
    memory->create(rgtall,maxchunk,6,"gyration/chunk:rgtall");
    array = rgtall;
  } else {
    memory->create(rg,maxchunk,"gyration/chunk:rg");
    memory->create(rgall,maxchunk,"gyration/chunk:rgall");
    vector = rgall;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeGyrationChunk::memory_usage()
{
  double bytes = (bigint) maxchunk * 2 * sizeof(double);
  bytes += (double) maxchunk * 2*3 * sizeof(double);
  if (tensor) bytes += (double) maxchunk * 2*6 * sizeof(double);
  else bytes += (double) maxchunk * 2 * sizeof(double);
  return bytes;
}

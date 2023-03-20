// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_vcm_chunk.h"

#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "compute_chunk_atom.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{ONCE,NFREQ,EVERY};

/* ---------------------------------------------------------------------- */

ComputeVCMChunk::ComputeVCMChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idchunk(nullptr), massproc(nullptr), masstotal(nullptr), vcm(nullptr), vcmall(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal compute vcm/chunk command");

  array_flag = 1;
  size_array_cols = 3;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  // ID of compute chunk/atom

  idchunk = utils::strdup(arg[3]);

  ComputeVCMChunk::init();

  // chunk-based data

  nchunk = 1;
  maxchunk = 0;
  allocate();

  firstflag = massneed = 1;
}

/* ---------------------------------------------------------------------- */

ComputeVCMChunk::~ComputeVCMChunk()
{
  delete [] idchunk;
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);
}

/* ---------------------------------------------------------------------- */

void ComputeVCMChunk::init()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute < 0)
    error->all(FLERR,"Chunk/atom compute does not exist for compute vcm/chunk");
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->compute[icompute]);
  if (strcmp(cchunk->style,"chunk/atom") != 0)
    error->all(FLERR,"Compute vcm/chunk does not use chunk/atom compute");
}

/* ---------------------------------------------------------------------- */

void ComputeVCMChunk::setup()
{
  // one-time calculation of per-chunk mass
  // done in setup, so that ComputeChunkAtom::setup() is already called

  if (firstflag && cchunk->idsflag == ONCE) {
    compute_array();
    firstflag = massneed = 0;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeVCMChunk::compute_array()
{
  int index;
  double massone;

  invoked_array = update->ntimestep;

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;

  if (nchunk > maxchunk) allocate();
  size_array_rows = nchunk;

  // zero local per-chunk values

  for (int i = 0; i < nchunk; i++)
    vcm[i][0] = vcm[i][1] = vcm[i][2] = 0.0;
  if (massneed)
    for (int i = 0; i < nchunk; i++) massproc[i] = 0.0;

  // compute VCM for each chunk

  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i]-1;
      if (index < 0) continue;
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      vcm[index][0] += v[i][0] * massone;
      vcm[index][1] += v[i][1] * massone;
      vcm[index][2] += v[i][2] * massone;
      if (massneed) massproc[index] += massone;
    }

  MPI_Allreduce(&vcm[0][0],&vcmall[0][0],3*nchunk,MPI_DOUBLE,MPI_SUM,world);
  if (massneed)
    MPI_Allreduce(massproc,masstotal,nchunk,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      vcmall[i][0] /= masstotal[i];
      vcmall[i][1] /= masstotal[i];
      vcmall[i][2] /= masstotal[i];
    } else vcmall[i][0] = vcmall[i][1] = vcmall[i][2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   lock methods: called by fix ave/time
   these methods ensure vector/array size is locked for Nfreq epoch
     by passing lock info along to compute chunk/atom
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   increment lock counter
------------------------------------------------------------------------- */

void ComputeVCMChunk::lock_enable()
{
  cchunk->lockcount++;
}

/* ----------------------------------------------------------------------
   decrement lock counter in compute chunk/atom, it if still exists
------------------------------------------------------------------------- */

void ComputeVCMChunk::lock_disable()
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

int ComputeVCMChunk::lock_length()
{
  nchunk = cchunk->setup_chunks();
  return nchunk;
}

/* ----------------------------------------------------------------------
   set the lock from startstep to stopstep
------------------------------------------------------------------------- */

void ComputeVCMChunk::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  cchunk->lock(fixptr,startstep,stopstep);
}

/* ----------------------------------------------------------------------
   unset the lock
------------------------------------------------------------------------- */

void ComputeVCMChunk::unlock(Fix *fixptr)
{
  cchunk->unlock(fixptr);
}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeVCMChunk::allocate()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);
  maxchunk = nchunk;
  memory->create(massproc,maxchunk,"vcm/chunk:massproc");
  memory->create(masstotal,maxchunk,"vcm/chunk:masstotal");
  memory->create(vcm,maxchunk,3,"vcm/chunk:vcm");
  memory->create(vcmall,maxchunk,3,"vcm/chunk:vcmall");
  array = vcmall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeVCMChunk::memory_usage()
{
  double bytes = (bigint) maxchunk * 2 * sizeof(double);
  bytes += (double) maxchunk * 2*3 * sizeof(double);
  return bytes;
}

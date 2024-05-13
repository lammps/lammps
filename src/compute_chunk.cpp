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

#include "compute_chunk.h"

#include "compute_chunk_atom.h"
#include "error.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;
/* ---------------------------------------------------------------------- */

ComputeChunk::ComputeChunk(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), idchunk(nullptr), cchunk(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, std::string("compute ") + style, error);

  // ID of compute chunk/atom

  idchunk = utils::strdup(arg[3]);

  ComputeChunk::init();

  // chunk-based data

  nchunk = 1;
  maxchunk = 0;
  firstflag = 1;
  massneed = 1;
}

/* ---------------------------------------------------------------------- */

ComputeChunk::~ComputeChunk()
{
  delete[] idchunk;
}

/* ---------------------------------------------------------------------- */

void ComputeChunk::init()
{
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (!cchunk)
    error->all(FLERR, "Chunk/atom compute {} does not exist or is incorrect style for compute {}",
               idchunk, style);
}

/* ----------------------------------------------------------------------
   common code used by all /chunk computes
------------------------------------------------------------------------- */

void ComputeChunk::compute_vector()
{
  invoked_vector = update->ntimestep;

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();

  if (nchunk > maxchunk) allocate();
  size_vector = nchunk;
}

/* ---------------------------------------------------------------------- */

void ComputeChunk::compute_array()
{
  invoked_array = update->ntimestep;

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();

  if (nchunk > maxchunk) allocate();
  size_array_rows = nchunk;
}

/* ----------------------------------------------------------------------
   lock methods: called by fix ave/time
   these methods ensure vector/array size is locked for Nfreq epoch
     by passing lock info along to compute chunk/atom
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   increment lock counter
------------------------------------------------------------------------- */

void ComputeChunk::lock_enable()
{
  cchunk->lockcount++;
}

/* ----------------------------------------------------------------------
   decrement lock counter in compute chunk/atom, it if still exists
------------------------------------------------------------------------- */

void ComputeChunk::lock_disable()
{
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (cchunk) cchunk->lockcount--;
}

/* ----------------------------------------------------------------------
   calculate and return # of chunks = length of vector/array
------------------------------------------------------------------------- */

int ComputeChunk::lock_length()
{
  nchunk = cchunk->setup_chunks();
  return nchunk;
}

/* ----------------------------------------------------------------------
   set the lock from startstep to stopstep
------------------------------------------------------------------------- */

void ComputeChunk::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  cchunk->lock(fixptr, startstep, stopstep);
}

/* ----------------------------------------------------------------------
   unset the lock
------------------------------------------------------------------------- */

void ComputeChunk::unlock(Fix *fixptr)
{
  cchunk->unlock(fixptr);
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeChunk::memory_usage()
{
  return sizeof(ComputeChunk) + sizeof(ComputeChunkAtom);
}

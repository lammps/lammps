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

#include "compute_gyration_chunk.h"

#include "atom.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGyrationChunk::ComputeGyrationChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), massproc(nullptr), masstotal(nullptr), com(nullptr),
    comall(nullptr), rg(nullptr), rgall(nullptr), rgt(nullptr), rgtall(nullptr)
{
  ComputeGyrationChunk::init();

  // optional args

  tensor = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "tensor") == 0) {
      tensor = 1;
      iarg++;
    } else
      error->all(FLERR, "Illegal compute gyration/chunk command");
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

  ComputeGyrationChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeGyrationChunk::~ComputeGyrationChunk()
{
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

void ComputeGyrationChunk::compute_vector()
{
  int i, index;
  double dx, dy, dz, massone;
  double unwrap[3];

  ComputeChunk::compute_vector();

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
      index = ichunk[i] - 1;
      if (index < 0) continue;
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      rg[index] += (dx * dx + dy * dy + dz * dz) * massone;
    }

  MPI_Allreduce(rg, rgall, nchunk, MPI_DOUBLE, MPI_SUM, world);

  for (i = 0; i < nchunk; i++)
    if (masstotal[i] > 0.0) rgall[i] = sqrt(rgall[i] / masstotal[i]);
}

/* ---------------------------------------------------------------------- */

void ComputeGyrationChunk::compute_array()
{
  int i, j, index;
  double dx, dy, dz, massone;
  double unwrap[3];

  ComputeChunk::compute_array();

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
      index = ichunk[i] - 1;
      if (index < 0) continue;
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      rgt[index][0] += dx * dx * massone;
      rgt[index][1] += dy * dy * massone;
      rgt[index][2] += dz * dz * massone;
      rgt[index][3] += dx * dy * massone;
      rgt[index][4] += dx * dz * massone;
      rgt[index][5] += dy * dz * massone;
    }

  if (nchunk) MPI_Allreduce(&rgt[0][0], &rgtall[0][0], nchunk * 6, MPI_DOUBLE, MPI_SUM, world);

  for (i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      for (j = 0; j < 6; j++) rgtall[i][j] = rgtall[i][j] / masstotal[i];
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

  int *ichunk = cchunk->ichunk;

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
      index = ichunk[i] - 1;
      if (index < 0) continue;
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      domain->unmap(x[i], image[i], unwrap);
      massproc[index] += massone;
      com[index][0] += unwrap[0] * massone;
      com[index][1] += unwrap[1] * massone;
      com[index][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(massproc, masstotal, nchunk, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&com[0][0], &comall[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);

  for (int i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      comall[i][0] /= masstotal[i];
      comall[i][1] /= masstotal[i];
      comall[i][2] /= masstotal[i];
    }
  }
}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeGyrationChunk::allocate()
{
  ComputeChunk::allocate();
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(rg);
  memory->destroy(rgall);
  memory->destroy(rgt);
  memory->destroy(rgtall);
  maxchunk = nchunk;
  memory->create(massproc, maxchunk, "gyration/chunk:massproc");
  memory->create(masstotal, maxchunk, "gyration/chunk:masstotal");
  memory->create(com, maxchunk, 3, "gyration/chunk:com");
  memory->create(comall, maxchunk, 3, "gyration/chunk:comall");
  if (tensor) {
    memory->create(rgt, maxchunk, 6, "gyration/chunk:rgt");
    memory->create(rgtall, maxchunk, 6, "gyration/chunk:rgtall");
    array = rgtall;
  } else {
    memory->create(rg, maxchunk, "gyration/chunk:rg");
    memory->create(rgall, maxchunk, "gyration/chunk:rgall");
    vector = rgall;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeGyrationChunk::memory_usage()
{
  double bytes = ComputeChunk::memory_usage();
  bytes += (bigint) maxchunk * 2 * sizeof(double);
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  if (tensor)
    bytes += (double) maxchunk * 2 * 6 * sizeof(double);
  else
    bytes += (double) maxchunk * 2 * sizeof(double);
  return bytes;
}

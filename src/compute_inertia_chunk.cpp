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

#include "compute_inertia_chunk.h"

#include "atom.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeInertiaChunk::ComputeInertiaChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), massproc(nullptr), masstotal(nullptr), com(nullptr),
    comall(nullptr), inertia(nullptr), inertiaall(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute inertia/chunk command");

  array_flag = 1;
  size_array_cols = 6;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  ComputeInertiaChunk::init();
  ComputeInertiaChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeInertiaChunk::~ComputeInertiaChunk()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(inertia);
  memory->destroy(inertiaall);
}

/* ---------------------------------------------------------------------- */

void ComputeInertiaChunk::compute_array()
{
  int i, j, index;
  double dx, dy, dz, massone;
  double unwrap[3];

  ComputeChunk::compute_array();
  int *ichunk = cchunk->ichunk;

  // zero local per-chunk values

  for (i = 0; i < nchunk; i++) {
    massproc[i] = 0.0;
    com[i][0] = com[i][1] = com[i][2] = 0.0;
    for (j = 0; j < 6; j++) inertia[i][j] = 0.0;
  }

  // compute COM for each chunk

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

  for (i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      comall[i][0] /= masstotal[i];
      comall[i][1] /= masstotal[i];
      comall[i][2] /= masstotal[i];
    }
  }

  // compute inertia tensor for each chunk

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];
      inertia[index][0] += massone * (dy * dy + dz * dz);
      inertia[index][1] += massone * (dx * dx + dz * dz);
      inertia[index][2] += massone * (dx * dx + dy * dy);
      inertia[index][3] -= massone * dx * dy;
      inertia[index][4] -= massone * dy * dz;
      inertia[index][5] -= massone * dx * dz;
    }

  MPI_Allreduce(&inertia[0][0], &inertiaall[0][0], 6 * nchunk, MPI_DOUBLE, MPI_SUM, world);
}
/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeInertiaChunk::allocate()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(inertia);
  memory->destroy(inertiaall);
  maxchunk = nchunk;
  memory->create(massproc, maxchunk, "inertia/chunk:massproc");
  memory->create(masstotal, maxchunk, "inertia/chunk:masstotal");
  memory->create(com, maxchunk, 3, "inertia/chunk:com");
  memory->create(comall, maxchunk, 3, "inertia/chunk:comall");
  memory->create(inertia, maxchunk, 6, "inertia/chunk:inertia");
  memory->create(inertiaall, maxchunk, 6, "inertia/chunk:inertiaall");
  array = inertiaall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeInertiaChunk::memory_usage()
{
  double bytes = ComputeChunk::memory_usage();
  bytes += (bigint) maxchunk * 2 * sizeof(double);
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  bytes += (double) maxchunk * 2 * 6 * sizeof(double);
  return bytes;
}

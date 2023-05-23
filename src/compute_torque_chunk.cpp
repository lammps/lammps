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

#include "compute_torque_chunk.h"

#include "atom.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTorqueChunk::ComputeTorqueChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), massproc(nullptr), masstotal(nullptr), com(nullptr),
    comall(nullptr), torque(nullptr), torqueall(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute torque/chunk command");

  array_flag = 1;
  size_array_cols = 3;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  ComputeTorqueChunk::init();
  ComputeTorqueChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeTorqueChunk::~ComputeTorqueChunk()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(torque);
  memory->destroy(torqueall);
}

/* ---------------------------------------------------------------------- */

void ComputeTorqueChunk::compute_array()
{
  int i, index;
  double dx, dy, dz, massone;
  double unwrap[3];

  ComputeChunk::compute_array();
  int *ichunk = cchunk->ichunk;

  // zero local per-chunk values

  for (i = 0; i < nchunk; i++) {
    massproc[i] = 0.0;
    com[i][0] = com[i][1] = com[i][2] = 0.0;
    torque[i][0] = torque[i][1] = torque[i][2] = 0.0;
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

  // compute torque on each chunk

  double **f = atom->f;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];
      torque[index][0] += dy * f[i][2] - dz * f[i][1];
      torque[index][1] += dz * f[i][0] - dx * f[i][2];
      torque[index][2] += dx * f[i][1] - dy * f[i][0];
    }

  MPI_Allreduce(&torque[0][0], &torqueall[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);
}
/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeTorqueChunk::allocate()
{
  ComputeChunk::allocate();
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(torque);
  memory->destroy(torqueall);
  maxchunk = nchunk;
  memory->create(massproc, maxchunk, "torque/chunk:massproc");
  memory->create(masstotal, maxchunk, "torque/chunk:masstotal");
  memory->create(com, maxchunk, 3, "torque/chunk:com");
  memory->create(comall, maxchunk, 3, "torque/chunk:comall");
  memory->create(torque, maxchunk, 3, "torque/chunk:torque");
  memory->create(torqueall, maxchunk, 3, "torque/chunk:torqueall");
  array = torqueall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeTorqueChunk::memory_usage()
{
  double bytes = (double) maxchunk * 2 * sizeof(double) + ComputeChunk::memory_usage();
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  return bytes;
}

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

#include "atom.h"
#include "compute_chunk_atom.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

enum { ONCE, NFREQ, EVERY };

/* ---------------------------------------------------------------------- */

ComputeVCMChunk::ComputeVCMChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), massproc(nullptr), masstotal(nullptr), vcm(nullptr),
    vcmall(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute vcm/chunk command");

  array_flag = 1;
  size_array_cols = 3;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  ComputeVCMChunk::init();
  ComputeVCMChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeVCMChunk::~ComputeVCMChunk()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);
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

  ComputeChunk::compute_array();
  int *ichunk = cchunk->ichunk;

  // zero local per-chunk values

  for (int i = 0; i < nchunk; i++) vcm[i][0] = vcm[i][1] = vcm[i][2] = 0.0;
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
      index = ichunk[i] - 1;
      if (index < 0) continue;
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      vcm[index][0] += v[i][0] * massone;
      vcm[index][1] += v[i][1] * massone;
      vcm[index][2] += v[i][2] * massone;
      if (massneed) massproc[index] += massone;
    }

  MPI_Allreduce(&vcm[0][0], &vcmall[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);
  if (massneed) MPI_Allreduce(massproc, masstotal, nchunk, MPI_DOUBLE, MPI_SUM, world);

  for (int i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      vcmall[i][0] /= masstotal[i];
      vcmall[i][1] /= masstotal[i];
      vcmall[i][2] /= masstotal[i];
    } else
      vcmall[i][0] = vcmall[i][1] = vcmall[i][2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeVCMChunk::allocate()
{
  ComputeChunk::allocate();
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);
  maxchunk = nchunk;
  memory->create(massproc, maxchunk, "vcm/chunk:massproc");
  memory->create(masstotal, maxchunk, "vcm/chunk:masstotal");
  memory->create(vcm, maxchunk, 3, "vcm/chunk:vcm");
  memory->create(vcmall, maxchunk, 3, "vcm/chunk:vcmall");
  array = vcmall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeVCMChunk::memory_usage()
{
  double bytes = (double) maxchunk * 2 * sizeof(double) + ComputeChunk::memory_usage();
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  return bytes;
}

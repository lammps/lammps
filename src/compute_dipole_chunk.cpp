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

#include "compute_dipole_chunk.h"

#include "atom.h"
#include "comm.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_special.h"
#include "memory.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathSpecial;

enum { MASSCENTER, GEOMCENTER };

/* ---------------------------------------------------------------------- */

ComputeDipoleChunk::ComputeDipoleChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), massproc(nullptr), masstotal(nullptr), chrgproc(nullptr),
    chrgtotal(nullptr), com(nullptr), comall(nullptr), dipole(nullptr), dipoleall(nullptr)
{
  if ((narg != 4) && (narg != 5)) error->all(FLERR, "Illegal compute dipole/chunk command");

  array_flag = 1;
  size_array_cols = 4;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  usecenter = MASSCENTER;

  if (narg == 5) {
    if (strncmp(arg[4], "geom", 4) == 0)
      usecenter = GEOMCENTER;
    else if (strcmp(arg[4], "mass") == 0)
      usecenter = MASSCENTER;
    else
      error->all(FLERR, "Illegal compute dipole/chunk command");
  }

  ComputeDipoleChunk::init();
  ComputeDipoleChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeDipoleChunk::~ComputeDipoleChunk()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(chrgproc);
  memory->destroy(chrgtotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(dipole);
  memory->destroy(dipoleall);
}

/* ---------------------------------------------------------------------- */

void ComputeDipoleChunk::init()
{
  ComputeChunk::init();

  if ((force->pair_match("tip4p/", 0) != nullptr) && (comm->me == 0))
    error->warning(FLERR, "Dipole moments may be incorrect when sing a TIP4P pair style");
}

/* ---------------------------------------------------------------------- */

void ComputeDipoleChunk::compute_array()
{
  int i, index;
  double massone;
  double unwrap[3];

  ComputeChunk::compute_array();
  int *ichunk = cchunk->ichunk;

  // zero local per-chunk values

  for (i = 0; i < nchunk; i++) {
    massproc[i] = chrgproc[i] = 0.0;
    com[i][0] = com[i][1] = com[i][2] = 0.0;
    dipole[i][0] = dipole[i][1] = dipole[i][2] = dipole[i][3] = 0.0;
  }

  // compute COM for each chunk

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *q = atom->q;
  double **mu = atom->mu;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      if (usecenter == MASSCENTER) {
        if (rmass)
          massone = rmass[i];
        else
          massone = mass[type[i]];
      } else
        massone = 1.0;    // usecenter == GEOMCENTER

      domain->unmap(x[i], image[i], unwrap);
      massproc[index] += massone;
      if (atom->q_flag) chrgproc[index] += atom->q[i];
      com[index][0] += unwrap[0] * massone;
      com[index][1] += unwrap[1] * massone;
      com[index][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(massproc, masstotal, nchunk, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(chrgproc, chrgtotal, nchunk, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&com[0][0], &comall[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);

  for (i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      comall[i][0] /= masstotal[i];
      comall[i][1] /= masstotal[i];
      comall[i][2] /= masstotal[i];
    }
  }

  // compute dipole for each chunk

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      domain->unmap(x[i], image[i], unwrap);
      if (atom->q_flag) {
        dipole[index][0] += q[i] * unwrap[0];
        dipole[index][1] += q[i] * unwrap[1];
        dipole[index][2] += q[i] * unwrap[2];
      }
      if (atom->mu_flag) {
        dipole[index][0] += mu[i][0];
        dipole[index][1] += mu[i][1];
        dipole[index][2] += mu[i][2];
      }
    }
  }

  MPI_Allreduce(&dipole[0][0], &dipoleall[0][0], 4 * nchunk, MPI_DOUBLE, MPI_SUM, world);

  for (i = 0; i < nchunk; i++) {
    // correct for position dependence with charged chunks
    dipoleall[i][0] -= chrgtotal[i] * comall[i][0];
    dipoleall[i][1] -= chrgtotal[i] * comall[i][1];
    dipoleall[i][2] -= chrgtotal[i] * comall[i][2];
    // compute total dipole moment
    dipoleall[i][3] =
        sqrt(square(dipoleall[i][0]) + square(dipoleall[i][1]) + square(dipoleall[i][2]));
  }
}
/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeDipoleChunk::allocate()
{
  ComputeChunk::allocate();
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(chrgproc);
  memory->destroy(chrgtotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(dipole);
  memory->destroy(dipoleall);
  maxchunk = nchunk;
  memory->create(massproc, maxchunk, "dipole/chunk:massproc");
  memory->create(masstotal, maxchunk, "dipole/chunk:masstotal");
  memory->create(chrgproc, maxchunk, "dipole/chunk:chrgproc");
  memory->create(chrgtotal, maxchunk, "dipole/chunk:chrgtotal");
  memory->create(com, maxchunk, 3, "dipole/chunk:com");
  memory->create(comall, maxchunk, 3, "dipole/chunk:comall");
  memory->create(dipole, maxchunk, 4, "dipole/chunk:dipole");
  memory->create(dipoleall, maxchunk, 4, "dipole/chunk:dipoleall");
  array = dipoleall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeDipoleChunk::memory_usage()
{
  double bytes = ComputeChunk::memory_usage();
  bytes += (bigint) maxchunk * 2 * sizeof(double);
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  bytes += (double) maxchunk * 2 * 4 * sizeof(double);
  return bytes;
}

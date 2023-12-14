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

#include "compute_msd_chunk.h"

#include "atom.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "fix_store_global.h"
#include "group.h"
#include "memory.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMSDChunk::ComputeMSDChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), id_fix(nullptr), fix(nullptr), massproc(nullptr),
    masstotal(nullptr), com(nullptr), comall(nullptr), msd(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute msd/chunk command");

  array_flag = 1;
  size_array_cols = 4;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  ComputeMSDChunk::init();

  // create a new fix STORE style for reference positions
  // id = compute-ID + COMPUTE_STORE, fix group = compute group
  // do not know size of array at this point, just allocate 1x1 array
  // fix creation must be done now so that a restart run can
  //   potentially re-populate the fix array (and change it to correct size)
  // otherwise size reset and init will be done in setup()

  id_fix = utils::strdup(std::string(id) + "_COMPUTE_STORE");
  fix = dynamic_cast<FixStoreGlobal *>(
      modify->add_fix(fmt::format("{} {} STORE/GLOBAL 1 1", id_fix, group->names[igroup])));
}

/* ---------------------------------------------------------------------- */

ComputeMSDChunk::~ComputeMSDChunk()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete[] id_fix;
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(msd);
}

/* ---------------------------------------------------------------------- */

void ComputeMSDChunk::init()
{
  ComputeChunk::init();

  // set fix which stores reference atom coords
  // if firstflag, will be created in setup()

  if (!firstflag) {
    fix = dynamic_cast<FixStoreGlobal *>(modify->get_fix_by_id(id_fix));
    if (!fix) error->all(FLERR, "Could not find compute msd/chunk fix with ID {}", id_fix);
  }
}

/* ----------------------------------------------------------------------
   compute initial COM for each chunk
   only once on timestep compute is defined, when firstflag = 1
------------------------------------------------------------------------- */

void ComputeMSDChunk::setup()
{
  if (!firstflag) return;
  compute_array();
  firstflag = 0;

  // if fix->astore is already correct size, restart file set it up
  // otherwise reset its size now and initialize to current COM

  if (fix->nrow == nchunk && fix->ncol == 3) return;
  fix->reset_global(nchunk, 3);

  double **cominit = fix->astore;
  for (int i = 0; i < nchunk; i++) {
    cominit[i][0] = comall[i][0];
    cominit[i][1] = comall[i][1];
    cominit[i][2] = comall[i][2];
    msd[i][0] = msd[i][1] = msd[i][2] = msd[i][3] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeMSDChunk::compute_array()
{
  int index;
  double massone;
  double unwrap[3];

  int oldnchunk = nchunk;
  ComputeChunk::compute_array();
  int *ichunk = cchunk->ichunk;

  // first time call, allocate per-chunk arrays
  // thereafter, require nchunk remain the same

  if (!firstflag && (oldnchunk != nchunk))
    error->all(FLERR, "Compute msd/chunk nchunk is not static");

  // zero local per-chunk values

  for (int i = 0; i < nchunk; i++) {
    massproc[i] = 0.0;
    com[i][0] = com[i][1] = com[i][2] = 0.0;
  }

  // compute current COM for each chunk

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

  // MSD is difference between current and initial COM
  // cominit is initilialized by setup() when firstflag is set

  if (firstflag) return;

  double dx, dy, dz;
  double **cominit = fix->astore;

  for (int i = 0; i < nchunk; i++) {
    dx = comall[i][0] - cominit[i][0];
    dy = comall[i][1] - cominit[i][1];
    dz = comall[i][2] - cominit[i][2];
    msd[i][0] = dx * dx;
    msd[i][1] = dy * dy;
    msd[i][2] = dz * dz;
    msd[i][3] = dx * dx + dy * dy + dz * dz;
  }
}

/* ----------------------------------------------------------------------
   one-time allocate of per-chunk arrays
------------------------------------------------------------------------- */

void ComputeMSDChunk::allocate()
{
  ComputeChunk::allocate();
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(msd);

  memory->create(massproc, nchunk, "msd/chunk:massproc");
  memory->create(masstotal, nchunk, "msd/chunk:masstotal");
  memory->create(com, nchunk, 3, "msd/chunk:com");
  memory->create(comall, nchunk, 3, "msd/chunk:comall");
  memory->create(msd, nchunk, 4, "msd/chunk:msd");
  array = msd;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeMSDChunk::memory_usage()
{
  double bytes = ComputeChunk::memory_usage();
  bytes += (bigint) nchunk * 2 * sizeof(double);
  bytes += (double) nchunk * 2 * 3 * sizeof(double);
  bytes += (double) nchunk * 4 * sizeof(double);
  return bytes;
}

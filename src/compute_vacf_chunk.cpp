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

#include "compute_vacf_chunk.h"

#include "atom.h"
#include "compute_chunk_atom.h"
#include "error.h"
#include "fix_store_global.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeVACFChunk::ComputeVACFChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), id_fix(nullptr), fix(nullptr), massproc(nullptr),
    masstotal(nullptr), vcm(nullptr), vcmall(nullptr), vacf(nullptr)
{
  if (narg != 4) error->all(FLERR, "Incorrect number of arguments for compute vacf/chunk");

  vacfnchunk = 0;
  array_flag = 1;
  size_array_cols = 4;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  ComputeVACFChunk::init();

  // create a new fix STORE style for reference velocities
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

ComputeVACFChunk::~ComputeVACFChunk()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete[] id_fix;
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);
  memory->destroy(vacf);
}

/* ---------------------------------------------------------------------- */

void ComputeVACFChunk::init()
{
  ComputeChunk::init();

  // set fix which stores reference atom coords
  // if firstflag, will be created in setup()

  if (!firstflag) {
    fix = dynamic_cast<FixStoreGlobal *>(modify->get_fix_by_id(id_fix));
    if (!fix) error->all(FLERR, "Could not find compute vacf/chunk fix with ID {}", id_fix);
  }
}

/* ----------------------------------------------------------------------
   compute initial VCM for each chunk
   only once on timestep compute is defined, when firstflag = 1
------------------------------------------------------------------------- */

void ComputeVACFChunk::setup()
{
  if (!firstflag) return;
  compute_array();
  firstflag = 0;

  // if fix->astore is already correct size, restart file set it up
  // otherwise reset its size now and initialize to current VCM

  if (fix->nrow == nchunk && fix->ncol == 3) return;
  fix->reset_global(nchunk, 3);

  double **vcminit = fix->astore;
  for (int i = 0; i < nchunk; i++) {
    vcminit[i][0] = vcmall[i][0];
    vcminit[i][1] = vcmall[i][1];
    vcminit[i][2] = vcmall[i][2];
    vacf[i][0] = vacf[i][1] = vacf[i][2] = vacf[i][3] = 1.0;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeVACFChunk::compute_array()
{
  invoked_array = update->ntimestep;

  int index;
  double massone;

  ComputeChunk::compute_array();
  int *ichunk = cchunk->ichunk;

  // first time call, allocate per-chunk arrays
  // thereafter, require nchunk remain the same

  if (firstflag)
    vacfnchunk = nchunk;
  else if (vacfnchunk != nchunk)
    error->all(FLERR, Error::NOLASTLINE, "Compute vacf/chunk nchunk is not static");

  // zero local per-chunk values

  for (int i = 0; i < nchunk; i++) {
    massproc[i] = 0.0;
    vcm[i][0] = vcm[i][1] = vcm[i][2] = 0.0;
  }

  // compute current VCM for each chunk

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
      massproc[index] += massone;
      vcm[index][0] += v[i][0] * massone;
      vcm[index][1] += v[i][1] * massone;
      vcm[index][2] += v[i][2] * massone;
    }

  MPI_Allreduce(massproc, masstotal, nchunk, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&vcm[0][0], &vcmall[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);

  for (int i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      vcmall[i][0] /= masstotal[i];
      vcmall[i][1] /= masstotal[i];
      vcmall[i][2] /= masstotal[i];
    }
  }

  // VACF is dot product between current and initial VCM
  // vcminit is initilialized by setup() when firstflag is set

  if (firstflag) return;

  double vxsq, vysq, vzsq;
  double **vcminit = fix->astore;

  for (int i = 0; i < nchunk; i++) {
    vxsq = vcmall[i][0] * vcminit[i][0];
    vysq = vcmall[i][1] * vcminit[i][1];
    vzsq = vcmall[i][2] * vcminit[i][2];
    vacf[i][0] = vxsq;
    vacf[i][1] = vysq;
    vacf[i][2] = vzsq;
    vacf[i][3] = vxsq + vysq + vzsq;
  }
}

/* ----------------------------------------------------------------------
   one-time allocate of per-chunk arrays
------------------------------------------------------------------------- */

void ComputeVACFChunk::allocate()
{
  ComputeChunk::allocate();
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);
  memory->destroy(vacf);

  memory->create(massproc, nchunk, "vacf/chunk:massproc");
  memory->create(masstotal, nchunk, "vacf/chunk:masstotal");
  memory->create(vcm, nchunk, 3, "vacf/chunk:vcm");
  memory->create(vcmall, nchunk, 3, "vacf/chunk:vcmall");
  memory->create(vacf, nchunk, 4, "vacf/chunk:vacf");
  array = vacf;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeVACFChunk::memory_usage()
{
  double bytes = ComputeChunk::memory_usage();
  bytes += (bigint) nchunk * 2 * sizeof(double);
  bytes += (double) nchunk * 2 * 3 * sizeof(double);
  bytes += (double) nchunk * 4 * sizeof(double);
  return bytes;
}

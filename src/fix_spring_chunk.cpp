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

#include "fix_spring_chunk.h"

#include "atom.h"
#include "comm.h"
#include "compute_chunk_atom.h"
#include "compute_com_chunk.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1.0e-10

/* ---------------------------------------------------------------------- */

FixSpringChunk::FixSpringChunk(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  idchunk(nullptr), idcom(nullptr), com0(nullptr), fcom(nullptr)
{
  if (narg != 6) utils::missing_cmd_args(FLERR, "fix spring/chunk", error);

  restart_global = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  k_spring = utils::numeric(FLERR,arg[3],false,lmp);

  idchunk = utils::strdup(arg[4]);
  idcom = utils::strdup(arg[5]);

  esprings = 0.0;
  nchunk = 0;
}

/* ---------------------------------------------------------------------- */

FixSpringChunk::~FixSpringChunk()
{
  memory->destroy(com0);
  memory->destroy(fcom);

  // decrement lock counter in compute chunk/atom, it if still exists

  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (cchunk) {
    cchunk->unlock(this);
    cchunk->lockcount--;
  }

  delete[] idchunk;
  delete[] idcom;
}

/* ---------------------------------------------------------------------- */

int FixSpringChunk::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpringChunk::init()
{
  // current indices for idchunk and idcom

  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (!cchunk)
    error->all(FLERR,"Chunk/atom compute {} does not exist or is not chunk/atom style", idchunk);

  ccom = dynamic_cast<ComputeCOMChunk *>(modify->get_compute_by_id(idcom));
  if (!ccom)
    error->all(FLERR,"Com/chunk compute {} does not exist or is not com/chunk style", idcom);

  // check that idchunk is consistent with ccom->idchunk

  if (ccom && (strcmp(idchunk,ccom->idchunk) != 0))
    error->all(FLERR,"Fix spring/chunk chunk ID {} not the same as compute com/chunk chunk ID {}",
               idchunk, ccom->idchunk);

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringChunk::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringChunk::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringChunk::post_force(int /*vflag*/)
{
  int i,m;
  double dx,dy,dz,r;

  // check if first time cchunk will be queried via ccom
  // if so, lock idchunk for as long as this fix is in place
  // will be unlocked in destructor
  // necessary b/c this fix stores original COM

  if (com0 == nullptr) cchunk->lock(this,update->ntimestep,-1);

  // calculate current centers of mass for each chunk
  // extract pointers from idchunk and idcom

  ccom->compute_array();

  nchunk = cchunk->nchunk;
  int *ichunk = cchunk->ichunk;
  double *masstotal = ccom->masstotal;
  double **com = ccom->array;

  // check if first time cchunk was queried via ccom
  // if so, allocate com0,fcom and store initial COM

  if (com0 == nullptr) {
    memory->create(com0,nchunk,3,"spring/chunk:com0");
    memory->create(fcom,nchunk,3,"spring/chunk:fcom");

    for (m = 0; m < nchunk; m++) {
      com0[m][0] = com[m][0];
      com0[m][1] = com[m][1];
      com0[m][2] = com[m][2];
    }
  }

  // calculate fcom = force on each COM, divided by masstotal

  esprings = 0.0;
  for (m = 0; m < nchunk; m++) {
    dx = com[m][0] - com0[m][0];
    dy = com[m][1] - com0[m][1];
    dz = com[m][2] - com0[m][2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    r = MAX(r,SMALL);

    if (masstotal[m]) {
      fcom[m][0] = k_spring*dx/r / masstotal[m];
      fcom[m][1] = k_spring*dy/r / masstotal[m];
      fcom[m][2] = k_spring*dz/r / masstotal[m];
      esprings += 0.5*k_spring*r*r;
    }
  }

  // apply restoring force to atoms in each chunk

  double **f = atom->f;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double massone;

  if (rmass) {
    for (i = 0; i < nlocal; i++) {
      m = ichunk[i]-1;
      if (m < 0) continue;
      massone = rmass[i];
      f[i][0] -= fcom[m][0]*massone;
      f[i][1] -= fcom[m][1]*massone;
      f[i][2] -= fcom[m][2]*massone;
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      m = ichunk[i]-1;
      if (m < 0) continue;
      massone = mass[type[i]];
      f[i][0] -= fcom[m][0]*massone;
      f[i][1] -= fcom[m][1]*massone;
      f[i][2] -= fcom[m][2]*massone;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringChunk::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringChunk::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   writ number of chunks and position of original COM into restart
------------------------------------------------------------------------- */

void FixSpringChunk::write_restart(FILE *fp)
{
  double n = nchunk;

  if (comm->me == 0) {
    int size = (3*n+1) * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&n,sizeof(double),1,fp);
    fwrite(&com0[0][0],3*sizeof(double),nchunk,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixSpringChunk::restart(char *buf)
{
  auto list = (double *) buf;
  int n = list[0];

  memory->destroy(com0);
  memory->destroy(fcom);

  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(idchunk));
  if (!cchunk)
    error->all(FLERR,"Chunk/atom compute {} does not exist or is not chunk/atom style", idchunk);

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  memory->create(com0,nchunk,3,"spring/chunk:com0");
  memory->create(fcom,nchunk,3,"spring/chunk:fcom");

  if (n != nchunk) {
    if (comm->me == 0)
      error->warning(FLERR,"Number of chunks changed from {} to {}. Cannot use restart", n, nchunk);
    memory->destroy(com0);
    memory->destroy(fcom);
    nchunk = 1;
  } else {
    cchunk->lock(this,update->ntimestep,-1);
    memcpy(&com0[0][0],list+1,3*n*sizeof(double));
  }
}

/* ----------------------------------------------------------------------
   energy of stretched spring
------------------------------------------------------------------------- */

double FixSpringChunk::compute_scalar()
{
  return esprings;
}

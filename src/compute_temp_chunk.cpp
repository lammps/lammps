/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "compute_temp_chunk.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "modify.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempChunk::ComputeTempChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute temp/chunk command");

  vector_flag = 1;
  size_vector = 1;
  size_vector_variable = 1;
  extvector = 0;

  // ID of compute chunk/atom

  int n = strlen(arg[3]) + 1;
  idchunk = new char[n];
  strcpy(idchunk,arg[3]);

  biasflag = 0;
  init();

  // optional args

  comflag = 0;
  biasflag = 0;
  id_bias = NULL;
  adof = domain->dimension;
  cdof = 0.0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"com") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/chunk command");
      if (strcmp(arg[iarg+1],"yes") == 0) comflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) comflag = 0;
      else error->all(FLERR,"Illegal compute temp/chunk command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"bias") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/chunk command");
      biasflag = 1;
      int n = strlen(arg[iarg+1]) + 1;
      id_bias = new char[n];
      strcpy(id_bias,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"adof") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/chunk command");
      adof = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"cdof") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/chunk command");
      cdof = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute temp/chunk command");
  }

  if (comflag && biasflag)
    error->all(FLERR,"Cannot use both com and bias with compute temp/chunk");

  // error check on bias compute

  if (biasflag) {
    int i = modify->find_compute(id_bias);
    if (i < 0)
      error->all(FLERR,"Could not find compute ID for temperature bias");
    tbias = modify->compute[i];
    if (tbias->tempflag == 0)
      error->all(FLERR,"Bias compute does not calculate temperature");
    if (tbias->tempbias == 0)
      error->all(FLERR,"Bias compute does not calculate a velocity bias");
  }

  // chunk-based data

  nchunk = 1;
  maxchunk = 0;
  ke = keall = NULL;
  count = countall = NULL;
  massproc = masstotal = NULL;
  vcm = vcmall = NULL;
  allocate();
}

/* ---------------------------------------------------------------------- */

ComputeTempChunk::~ComputeTempChunk()
{
  delete [] idchunk;
  delete [] id_bias;
  memory->destroy(ke);
  memory->destroy(keall);
  memory->destroy(count);
  memory->destroy(countall);
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(vcm);
  memory->destroy(vcmall);
}

/* ---------------------------------------------------------------------- */

void ComputeTempChunk::init()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute < 0)
    error->all(FLERR,"Chunk/atom compute does not exist for "
               "compute temp/chunk");
  cchunk = (ComputeChunkAtom *) modify->compute[icompute];
  if (strcmp(cchunk->style,"chunk/atom") != 0)
    error->all(FLERR,"Compute temp/chunk does not use chunk/atom compute");

  if (biasflag) {
    int i = modify->find_compute(id_bias);
    if (i < 0)
      error->all(FLERR,"Could not find compute ID for temperature bias");
    tbias = modify->compute[i];
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTempChunk::compute_vector()
{
  int index;

  invoked_vector = update->ntimestep;

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;

  if (nchunk > maxchunk) allocate();

  // calculate COM velocity for each chunk

  if (comflag) vcm_compute();

  // remove velocity bias

  if (biasflag) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  // zero local per-chunk values

  for (int i = 0; i < nchunk; i++) {
    count[i] = 0;
    ke[i] = 0.0;
  }

  // compute temperature for each chunk
  // option for removing COM velocity

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  if (!comflag) {
    if (rmass) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          ke[index] += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * 
            rmass[i];
          count[index]++;
        }
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          ke[index] += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
            mass[type[i]];
          count[index]++;
      }
    }

  } else {
    double vx,vy,vz;
    if (rmass) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          vx = v[i][0] - vcmall[index][0];
          vy = v[i][1] - vcmall[index][1];
          vz = v[i][2] - vcmall[index][2];
          ke[index] += (vx*vx + vy*vy + vz*vz) * rmass[i];
          count[index]++;
        }
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          vx = v[i][0] - vcmall[index][0];
          vy = v[i][1] - vcmall[index][1];
          vz = v[i][2] - vcmall[index][2];
          ke[index] += (vx*vx + vy*vy + vz*vz) * mass[type[i]];
          count[index]++;
      }
    }
  }

  MPI_Allreduce(ke,keall,nchunk,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(count,countall,nchunk,MPI_INT,MPI_SUM,world);

  // restore velocity bias

  if (biasflag) tbias->restore_bias_all();

  // normalize temperatures by per-chunk DOF

  double dof,tfactor;
  double mvv2e = force->mvv2e;
  double boltz = force->boltz;

  for (int i = 0; i < nchunk; i++) {
    dof = cdof + adof*countall[i];
    if (dof > 0.0) tfactor = mvv2e / (dof * boltz);
    else tfactor = 0.0;
    keall[i] *= tfactor;
  }
}

/* ----------------------------------------------------------------------
   calculate velocity of COM for each chunk
------------------------------------------------------------------------- */

void ComputeTempChunk::vcm_compute()
{
  int index;
  double massone;

  int *ichunk = cchunk->ichunk;

  for (int i = 0; i < nchunk; i++) {
    vcm[i][0] = vcm[i][1] = vcm[i][2] = 0.0;
    massproc[i] = 0.0;
  }

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
      massproc[index] += massone;
    }

  MPI_Allreduce(&vcm[0][0],&vcmall[0][0],3*nchunk,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(massproc,masstotal,nchunk,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0; i < nchunk; i++) {
    vcmall[i][0] /= masstotal[i];
    vcmall[i][1] /= masstotal[i];
    vcmall[i][2] /= masstotal[i];
  }
}

/* ----------------------------------------------------------------------
   lock methods: called by fix ave/time
   these methods insure vector/array size is locked for Nfreq epoch
     by passing lock info along to compute chunk/atom
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   increment lock counter
------------------------------------------------------------------------- */

void ComputeTempChunk::lock_enable()
{
  cchunk->lockcount++;
}

/* ----------------------------------------------------------------------
   decrement lock counter in compute chunk/atom, it if still exists
------------------------------------------------------------------------- */

void ComputeTempChunk::lock_disable()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute >= 0) {
    cchunk = (ComputeChunkAtom *) modify->compute[icompute];
    cchunk->lockcount--;
  }
}

/* ----------------------------------------------------------------------
   calculate and return # of chunks = length of vector/array
------------------------------------------------------------------------- */

int ComputeTempChunk::lock_length()
{
  nchunk = cchunk->setup_chunks();
  return nchunk;
}

/* ----------------------------------------------------------------------
   set the lock from startstep to stopstep
------------------------------------------------------------------------- */

void ComputeTempChunk::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  cchunk->lock(fixptr,startstep,stopstep);
}

/* ----------------------------------------------------------------------
   unset the lock
------------------------------------------------------------------------- */

void ComputeTempChunk::unlock(Fix *fixptr)
{
  cchunk->unlock(fixptr);
}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeTempChunk::allocate()
{
  memory->destroy(ke);
  memory->destroy(keall);
  memory->destroy(count);
  memory->destroy(countall);
  size_vector = maxchunk = nchunk;
  memory->create(ke,maxchunk,"temp/chunk:ke");
  memory->create(keall,maxchunk,"temp/chunk:keall");
  memory->create(count,maxchunk,"temp/chunk:count");
  memory->create(countall,maxchunk,"temp/chunk:countall");
  vector = keall;

  if (comflag) {
    memory->destroy(massproc);
    memory->destroy(masstotal);
    memory->destroy(vcm);
    memory->destroy(vcmall);
    memory->create(massproc,maxchunk,"vcm/chunk:massproc");
    memory->create(masstotal,maxchunk,"vcm/chunk:masstotal");
    memory->create(vcm,maxchunk,3,"vcm/chunk:vcm");
    memory->create(vcmall,maxchunk,3,"vcm/chunk:vcmall");
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeTempChunk::memory_usage()
{
  double bytes = (bigint) maxchunk * 2 * sizeof(double);
  bytes = (bigint) maxchunk * 2 * sizeof(int);
  if (comflag) {
    bytes += (bigint) maxchunk * 2 * sizeof(double);
    bytes += (bigint) maxchunk * 2*3 * sizeof(double);
  }
  return bytes;
}

// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_temp_chunk.h"

#include "atom.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum{TEMP,KECOM,INTERNAL};

/* ---------------------------------------------------------------------- */

ComputeTempChunk::ComputeTempChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  which(nullptr), idchunk(nullptr), id_bias(nullptr), sum(nullptr), sumall(nullptr), count(nullptr),
  countall(nullptr), massproc(nullptr), masstotal(nullptr), vcm(nullptr), vcmall(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute temp/chunk command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  // ID of compute chunk/atom

  idchunk = utils::strdup(arg[3]);

  biasflag = 0;
  ComputeTempChunk::init();

  // optional per-chunk values

  nvalues = narg-4;
  which = new int[nvalues];
  nvalues = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"temp") == 0) which[nvalues] = TEMP;
    else if (strcmp(arg[iarg],"kecom") == 0) which[nvalues] = KECOM;
    else if (strcmp(arg[iarg],"internal") == 0) which[nvalues] = INTERNAL;
    else break;
    iarg++;
    nvalues++;
  }

  // optional args

  comflag = 0;
  biasflag = 0;
  id_bias = nullptr;
  adof = domain->dimension;
  cdof = 0.0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"com") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp/chunk command");
      comflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bias") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp/chunk command");
      biasflag = 1;
      id_bias = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"adof") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp/chunk command");
      adof = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"cdof") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp/chunk command");
      cdof = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute temp/chunk command");
  }

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

  // this compute only calculates a bias, if comflag is set
  // won't be two biases since comflag and biasflag cannot both be set

  if (comflag && biasflag)
    error->all(FLERR,"Cannot use both com and bias with compute temp/chunk");
  if (comflag) tempbias = 1;

  // vector data

  vector = new double[size_vector];

  // chunk-based data

  nchunk = 1;
  maxchunk = 0;

  if (nvalues)  {
    array_flag = 1;
    size_array_cols = nvalues;
    size_array_rows = 0;
    size_array_rows_variable = 1;
    extarray = 0;
  }

  allocate();
  comstep = -1;
}

/* ---------------------------------------------------------------------- */

ComputeTempChunk::~ComputeTempChunk()
{
  delete [] idchunk;
  delete [] which;
  delete [] id_bias;
  delete [] vector;
  memory->destroy(sum);
  memory->destroy(sumall);
  memory->destroy(count);
  memory->destroy(countall);
  memory->destroy(array);
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
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->compute[icompute]);
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

double ComputeTempChunk::compute_scalar()
{
  int i,index;

  invoked_scalar = update->ntimestep;

  // calculate chunk assignments,
  //   since only atoms in chunks contribute to global temperature
  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;

  if (nchunk > maxchunk) allocate();

  // remove velocity bias

  if (biasflag) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  // calculate COM velocity for each chunk
  // won't be invoked with bias also removed = 2 biases

  if (comflag && comstep != update->ntimestep) vcm_compute();

  // calculate global temperature, optionally removing COM velocity

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;
  int mycount = 0;

  if (!comflag) {
    if (rmass) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
            rmass[i];
          mycount++;
        }
    } else {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
            mass[type[i]];
          mycount++;
      }
    }

  } else {
    double vx,vy,vz;
    if (rmass) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          vx = v[i][0] - vcmall[index][0];
          vy = v[i][1] - vcmall[index][1];
          vz = v[i][2] - vcmall[index][2];
          t += (vx*vx + vy*vy + vz*vz) * rmass[i];
          mycount++;
        }
    } else {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          vx = v[i][0] - vcmall[index][0];
          vy = v[i][1] - vcmall[index][1];
          vz = v[i][2] - vcmall[index][2];
          t += (vx*vx + vy*vy + vz*vz) * mass[type[i]];
          mycount++;
      }
    }
  }

  // restore velocity bias

  if (biasflag) tbias->restore_bias_all();

  // final temperature

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  double rcount = mycount;
  double allcount;
  MPI_Allreduce(&rcount,&allcount,1,MPI_DOUBLE,MPI_SUM,world);

  double dof = nchunk*cdof + adof*allcount;
  double tfactor = 0.0;
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  if (dof < 0.0 && allcount > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempChunk::compute_vector()
{
  int i,index;

  invoked_vector = update->ntimestep;

  // calculate chunk assignments,
  //   since only atoms in chunks contribute to global temperature
  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;

  if (nchunk > maxchunk) allocate();

  // remove velocity bias

  if (biasflag) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  // calculate COM velocity for each chunk
  // won't be invoked with bias also removed = 2 biases

  if (comflag && comstep != update->ntimestep) vcm_compute();

  // calculate KE tensor, optionally removing COM velocity

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  if (!comflag) {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        index = ichunk[i]-1;
        if (index < 0) continue;
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        t[0] += massone * v[i][0]*v[i][0];
        t[1] += massone * v[i][1]*v[i][1];
        t[2] += massone * v[i][2]*v[i][2];
        t[3] += massone * v[i][0]*v[i][1];
        t[4] += massone * v[i][0]*v[i][2];
        t[5] += massone * v[i][1]*v[i][2];
      }
  } else {
    double vx,vy,vz;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        index = ichunk[i]-1;
        if (index < 0) continue;
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        vx = v[i][0] - vcmall[index][0];
        vy = v[i][1] - vcmall[index][1];
        vz = v[i][2] - vcmall[index][2];
        t[0] += massone * vx*vx;
        t[1] += massone * vy*vy;
        t[2] += massone * vz*vz;
        t[3] += massone * vx*vy;
        t[4] += massone * vx*vz;
        t[5] += massone * vy*vz;
      }
  }

  // restore velocity bias

  if (biasflag) tbias->restore_bias_all();

  // final KE

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ---------------------------------------------------------------------- */

void ComputeTempChunk::compute_array()
{
  invoked_array = update->ntimestep;

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();

  if (nchunk > maxchunk) allocate();
  size_array_rows = nchunk;

  // remove velocity bias

  if (biasflag) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  // calculate COM velocity for each chunk whether comflag set or not
  //   needed by some values even if comflag not set
  // important to do this after velocity bias is removed
  //   otherwise per-chunk values that use both v and vcm will be inconsistent

  if (comstep != update->ntimestep) vcm_compute();

  // compute each value

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == TEMP) temperature(i);
    else if (which[i] == KECOM) kecom(i);
    else if (which[i] == INTERNAL) internal(i);
  }

  // restore velocity bias

  if (biasflag) tbias->restore_bias_all();
}

/* ----------------------------------------------------------------------
   calculate velocity of COM for each chunk
------------------------------------------------------------------------- */

void ComputeTempChunk::vcm_compute()
{
  int i,index;
  double massone;

  // avoid re-computing VCM more than once per step

  comstep = update->ntimestep;

  int *ichunk = cchunk->ichunk;

  for (i = 0; i < nchunk; i++) {
    vcm[i][0] = vcm[i][1] = vcm[i][2] = 0.0;
    massproc[i] = 0.0;
  }

  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
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

  for (i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      vcmall[i][0] /= masstotal[i];
      vcmall[i][1] /= masstotal[i];
      vcmall[i][2] /= masstotal[i];
    } else {
      vcmall[i][0] = vcmall[i][1] = vcmall[i][2] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   temperature of each chunk
------------------------------------------------------------------------- */

void ComputeTempChunk::temperature(int icol)
{
  int i,index;
  int *ichunk = cchunk->ichunk;

  // zero local per-chunk values

  for (i = 0; i < nchunk; i++) {
    count[i] = 0;
    sum[i] = 0.0;
  }

  // per-chunk temperature, option for removing COM velocity

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  if (!comflag) {
    if (rmass) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          sum[index] += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
            rmass[i];
          count[index]++;
        }
    } else {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          sum[index] += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
            mass[type[i]];
          count[index]++;
      }
    }

  } else {
    double vx,vy,vz;
    if (rmass) {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          vx = v[i][0] - vcmall[index][0];
          vy = v[i][1] - vcmall[index][1];
          vz = v[i][2] - vcmall[index][2];
          sum[index] += (vx*vx + vy*vy + vz*vz) * rmass[i];
          count[index]++;
        }
    } else {
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          index = ichunk[i]-1;
          if (index < 0) continue;
          vx = v[i][0] - vcmall[index][0];
          vy = v[i][1] - vcmall[index][1];
          vz = v[i][2] - vcmall[index][2];
          sum[index] += (vx*vx + vy*vy + vz*vz) * mass[type[i]];
          count[index]++;
      }
    }
  }

  // sum across procs

  MPI_Allreduce(sum,sumall,nchunk,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(count,countall,nchunk,MPI_INT,MPI_SUM,world);

  // normalize temperatures by per-chunk DOF

  double dof,tfactor;
  double mvv2e = force->mvv2e;
  double boltz = force->boltz;

  for (i = 0; i < nchunk; i++) {
    dof = cdof + adof*countall[i];
    if (dof > 0.0) tfactor = mvv2e / (dof * boltz);
    else tfactor = 0.0;
    array[i][icol] = tfactor * sumall[i];
  }
}

/* ----------------------------------------------------------------------
   KE of entire chunk moving at VCM
------------------------------------------------------------------------- */

void ComputeTempChunk::kecom(int icol)
{
  int index;
  int *ichunk = cchunk->ichunk;

  // zero local per-chunk values

  for (int i = 0; i < nchunk; i++) sum[i] = 0.0;

  // per-chunk COM KE

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double vx,vy,vz;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        index = ichunk[i]-1;
        if (index < 0) continue;
        vx = vcmall[index][0];
        vy = vcmall[index][1];
        vz = vcmall[index][2];
        sum[index] += (vx*vx + vy*vy + vz*vz) * rmass[i];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        index = ichunk[i]-1;
        if (index < 0) continue;
        vx = vcmall[index][0];
        vy = vcmall[index][1];
        vz = vcmall[index][2];
        sum[index] += (vx*vx + vy*vy + vz*vz) * mass[type[i]];
      }
  }

  // sum across procs

  MPI_Allreduce(sum,sumall,nchunk,MPI_DOUBLE,MPI_SUM,world);

  double mvv2e = force->mvv2e;
  for (int i = 0; i < nchunk; i++)
    array[i][icol] = 0.5 * mvv2e * sumall[i];

}

/* ----------------------------------------------------------------------
   internal KE of each chunk around its VCM
   computed using per-atom velocities with chunk VCM subtracted off
------------------------------------------------------------------------- */

void ComputeTempChunk::internal(int icol)
{
  int index;
  int *ichunk = cchunk->ichunk;

  // zero local per-chunk values

  for (int i = 0; i < nchunk; i++) sum[i] = 0.0;

  // per-chunk internal KE

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double vx,vy,vz;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        index = ichunk[i]-1;
        if (index < 0) continue;
        vx = v[i][0] - vcmall[index][0];
        vy = v[i][1] - vcmall[index][1];
        vz = v[i][2] - vcmall[index][2];
        sum[index] += (vx*vx + vy*vy + vz*vz) * rmass[i];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        index = ichunk[i]-1;
        if (index < 0) continue;
        vx = v[i][0] - vcmall[index][0];
        vy = v[i][1] - vcmall[index][1];
        vz = v[i][2] - vcmall[index][2];
        sum[index] += (vx*vx + vy*vy + vz*vz) * mass[type[i]];
      }
  }

  // sum across procs

  MPI_Allreduce(sum,sumall,nchunk,MPI_DOUBLE,MPI_SUM,world);

  double mvv2e = force->mvv2e;
  for (int i = 0; i < nchunk; i++)
    array[i][icol] = 0.5 * mvv2e * sumall[i];
}

/* ----------------------------------------------------------------------
   bias methods: called by thermostats
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempChunk::remove_bias(int i, double *v)
{
  int index = cchunk->ichunk[i]-1;
  if (index < 0) return;
  v[0] -= vcmall[index][0];
  v[1] -= vcmall[index][1];
  v[2] -= vcmall[index][2];
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempChunk::remove_bias_all()
{
  int index;
  int *ichunk = cchunk->ichunk;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i]-1;
      if (index < 0) continue;
      v[i][0] -= vcmall[index][0];
      v[i][1] -= vcmall[index][1];
      v[i][2] -= vcmall[index][2];
    }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempChunk::restore_bias(int i, double *v)
{
  int index = cchunk->ichunk[i]-1;
  if (index < 0) return;
  v[0] += vcmall[index][0];
  v[1] += vcmall[index][1];
  v[2] += vcmall[index][2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempChunk::restore_bias_all()
{
  int index;
  int *ichunk = cchunk->ichunk;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i]-1;
      if (index < 0) continue;
      v[i][0] += vcmall[index][0];
      v[i][1] += vcmall[index][1];
      v[i][2] += vcmall[index][2];
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
    cchunk = dynamic_cast<ComputeChunkAtom *>(modify->compute[icompute]);
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
  memory->destroy(sum);
  memory->destroy(sumall);
  memory->destroy(count);
  memory->destroy(countall);
  memory->destroy(array);
  maxchunk = nchunk;
  memory->create(sum,maxchunk,"temp/chunk:sum");
  memory->create(sumall,maxchunk,"temp/chunk:sumall");
  memory->create(count,maxchunk,"temp/chunk:count");
  memory->create(countall,maxchunk,"temp/chunk:countall");
  memory->create(array,maxchunk,nvalues,"temp/chunk:array");

  if (comflag || nvalues) {
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
  bytes += (double) maxchunk * 2 * sizeof(int);
  bytes += (double) maxchunk * nvalues * sizeof(double);
  if (comflag || nvalues) {
    bytes += (double) maxchunk * 2 * sizeof(double);
    bytes += (double) maxchunk * 2*3 * sizeof(double);
  }
  return bytes;
}

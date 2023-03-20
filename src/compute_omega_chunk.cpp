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

#include "compute_omega_chunk.h"

#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "math_extra.h"
#include "math_eigen.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

ComputeOmegaChunk::ComputeOmegaChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idchunk(nullptr),massproc(nullptr),masstotal(nullptr),com(nullptr),comall(nullptr),
  inertia(nullptr),inertiaall(nullptr),angmom(nullptr),angmomall(nullptr),omega(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal compute omega/chunk command");

  array_flag = 1;
  size_array_cols = 3;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  // ID of compute chunk/atom

  idchunk = utils::strdup(arg[3]);

  ComputeOmegaChunk::init();

  // chunk-based data

  nchunk = 1;
  maxchunk = 0;
  allocate();
}

/* ---------------------------------------------------------------------- */

ComputeOmegaChunk::~ComputeOmegaChunk()
{
  delete [] idchunk;
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(angmom);
  memory->destroy(angmomall);
  memory->destroy(inertia);
  memory->destroy(inertiaall);
  memory->destroy(omega);
}

/* ---------------------------------------------------------------------- */

void ComputeOmegaChunk::init()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute < 0)
    error->all(FLERR,"Chunk/atom compute does not exist for "
               "compute omega/chunk");
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->compute[icompute]);
  if (strcmp(cchunk->style,"chunk/atom") != 0)
    error->all(FLERR,"Compute omega/chunk does not use chunk/atom compute");
}

/* ---------------------------------------------------------------------- */

void ComputeOmegaChunk::compute_array()
{
  int i,j,m,index;
  double dx,dy,dz,massone;
  double unwrap[3];

  invoked_array = update->ntimestep;

  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;

  if (nchunk > maxchunk) allocate();
  size_array_rows = nchunk;

  // zero local per-chunk values

  for (i = 0; i < nchunk; i++) {
    massproc[i] = 0.0;
    com[i][0] = com[i][1] = com[i][2] = 0.0;
    for (j = 0; j < 6; j++) inertia[i][j] = 0.0;
    angmom[i][0] = angmom[i][1] = angmom[i][2] = 0.0;
    omega[i][0] = omega[i][1] = omega[i][2] = 0.0;
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
      index = ichunk[i]-1;
      if (index < 0) continue;
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      domain->unmap(x[i],image[i],unwrap);
      massproc[index] += massone;
      com[index][0] += unwrap[0] * massone;
      com[index][1] += unwrap[1] * massone;
      com[index][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(massproc,masstotal,nchunk,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&com[0][0],&comall[0][0],3*nchunk,MPI_DOUBLE,MPI_SUM,world);

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
      index = ichunk[i]-1;
      if (index < 0) continue;
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];
      inertia[index][0] += massone * (dy*dy + dz*dz);
      inertia[index][1] += massone * (dx*dx + dz*dz);
      inertia[index][2] += massone * (dx*dx + dy*dy);
      inertia[index][3] -= massone * dx*dy;
      inertia[index][4] -= massone * dy*dz;
      inertia[index][5] -= massone * dx*dz;
    }

  MPI_Allreduce(&inertia[0][0],&inertiaall[0][0],6*nchunk,
                MPI_DOUBLE,MPI_SUM,world);

  // compute angmom for each chunk

  double **v = atom->v;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i]-1;
      if (index < 0) continue;
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      angmom[index][0] += massone * (dy*v[i][2] - dz*v[i][1]);
      angmom[index][1] += massone * (dz*v[i][0] - dx*v[i][2]);
      angmom[index][2] += massone * (dx*v[i][1] - dy*v[i][0]);
    }

  MPI_Allreduce(&angmom[0][0],&angmomall[0][0],3*nchunk,
                MPI_DOUBLE,MPI_SUM,world);

  // compute omega for each chunk

  double determinant,invdeterminant;
  double idiag[3],ex[3],ey[3],ez[3],cross[3];
  double ione[3][3],inverse[3][3],evectors[3][3];
  double *iall,*mall;

  for (m = 0; m < nchunk; m++) {

    // determinant = triple product of rows of inertia matrix

    iall = &inertiaall[m][0];
    determinant = iall[0] * (iall[1]*iall[2] - iall[4]*iall[4]) +
      iall[3] * (iall[4]*iall[5] - iall[3]*iall[2]) +
      iall[5] * (iall[3]*iall[4] - iall[1]*iall[5]);

    ione[0][0] = iall[0];
    ione[1][1] = iall[1];
    ione[2][2] = iall[2];
    ione[0][1] = ione[1][0] = iall[3];
    ione[1][2] = ione[2][1] = iall[4];
    ione[0][2] = ione[2][0] = iall[5];

    // non-singular I matrix
    // use L = Iw, inverting I to solve for w

    if (determinant > EPSILON) {
      inverse[0][0] = ione[1][1]*ione[2][2] - ione[1][2]*ione[2][1];
      inverse[0][1] = -(ione[0][1]*ione[2][2] - ione[0][2]*ione[2][1]);
      inverse[0][2] = ione[0][1]*ione[1][2] - ione[0][2]*ione[1][1];

      inverse[1][0] = -(ione[1][0]*ione[2][2] - ione[1][2]*ione[2][0]);
      inverse[1][1] = ione[0][0]*ione[2][2] - ione[0][2]*ione[2][0];
      inverse[1][2] = -(ione[0][0]*ione[1][2] - ione[0][2]*ione[1][0]);

      inverse[2][0] = ione[1][0]*ione[2][1] - ione[1][1]*ione[2][0];
      inverse[2][1] = -(ione[0][0]*ione[2][1] - ione[0][1]*ione[2][0]);
      inverse[2][2] = ione[0][0]*ione[1][1] - ione[0][1]*ione[1][0];

      invdeterminant = 1.0/determinant;
      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
          inverse[i][j] *= invdeterminant;

      mall = &angmomall[m][0];
      omega[m][0] = inverse[0][0]*mall[0] + inverse[0][1]*mall[1] +
        inverse[0][2]*mall[2];
      omega[m][1] = inverse[1][0]*mall[0] + inverse[1][1]*mall[1] +
        inverse[1][2]*mall[2];
      omega[m][2] = inverse[2][0]*mall[0] + inverse[2][1]*mall[1] +
        inverse[2][2]*mall[2];

    // handle each (nearly) singular I matrix
    // due to 2-atom chunk or linear molecule
    // use jacobi3() and angmom_to_omega() to calculate valid omega

    } else {
      int ierror = MathEigen::jacobi3(ione,idiag,evectors);
      if (ierror) error->all(FLERR,
                             "Insufficient Jacobi rotations for omega/chunk");

      ex[0] = evectors[0][0];
      ex[1] = evectors[1][0];
      ex[2] = evectors[2][0];
      ey[0] = evectors[0][1];
      ey[1] = evectors[1][1];
      ey[2] = evectors[2][1];
      ez[0] = evectors[0][2];
      ez[1] = evectors[1][2];
      ez[2] = evectors[2][2];

      // enforce 3 evectors as a right-handed coordinate system
      // flip 3rd vector if needed

      MathExtra::cross3(ex,ey,cross);
      if (MathExtra::dot3(cross,ez) < 0.0) MathExtra::negate3(ez);

      // if any principal moment < scaled EPSILON, set to 0.0

      double max;
      max = MAX(idiag[0],idiag[1]);
      max = MAX(max,idiag[2]);

      if (idiag[0] < EPSILON*max) idiag[0] = 0.0;
      if (idiag[1] < EPSILON*max) idiag[1] = 0.0;
      if (idiag[2] < EPSILON*max) idiag[2] = 0.0;

      // calculate omega using diagonalized inertia matrix

      MathExtra::angmom_to_omega(&angmomall[m][0],ex,ey,ez,idiag,&omega[m][0]);
    }
  }
}

/* ----------------------------------------------------------------------
   lock methods: called by fix ave/time
   these methods ensure vector/array size is locked for Nfreq epoch
     by passing lock info along to compute chunk/atom
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   increment lock counter
------------------------------------------------------------------------- */

void ComputeOmegaChunk::lock_enable()
{
  cchunk->lockcount++;
}

/* ----------------------------------------------------------------------
   decrement lock counter in compute chunk/atom, it if still exists
------------------------------------------------------------------------- */

void ComputeOmegaChunk::lock_disable()
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

int ComputeOmegaChunk::lock_length()
{
  nchunk = cchunk->setup_chunks();
  return nchunk;
}

/* ----------------------------------------------------------------------
   set the lock from startstep to stopstep
------------------------------------------------------------------------- */

void ComputeOmegaChunk::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  cchunk->lock(fixptr,startstep,stopstep);
}

/* ----------------------------------------------------------------------
   unset the lock
------------------------------------------------------------------------- */

void ComputeOmegaChunk::unlock(Fix *fixptr)
{
  cchunk->unlock(fixptr);
}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeOmegaChunk::allocate()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(inertia);
  memory->destroy(inertiaall);
  memory->destroy(angmom);
  memory->destroy(angmomall);
  memory->destroy(omega);
  maxchunk = nchunk;
  memory->create(massproc,maxchunk,"omega/chunk:massproc");
  memory->create(masstotal,maxchunk,"omega/chunk:masstotal");
  memory->create(com,maxchunk,3,"omega/chunk:com");
  memory->create(comall,maxchunk,3,"omega/chunk:comall");
  memory->create(inertia,maxchunk,6,"omega/chunk:inertia");
  memory->create(inertiaall,maxchunk,6,"omega/chunk:inertiaall");
  memory->create(angmom,maxchunk,3,"omega/chunk:angmom");
  memory->create(angmomall,maxchunk,3,"omega/chunk:angmomall");
  memory->create(omega,maxchunk,3,"omega/chunk:omega");
  array = omega;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeOmegaChunk::memory_usage()
{
  double bytes = (bigint) maxchunk * 2 * sizeof(double);
  bytes += (double) maxchunk * 2*3 * sizeof(double);
  bytes += (double) maxchunk * 2*6 * sizeof(double);
  bytes += (double) maxchunk * 2*3 * sizeof(double);
  bytes += (double) maxchunk * 3 * sizeof(double);
  return bytes;
}

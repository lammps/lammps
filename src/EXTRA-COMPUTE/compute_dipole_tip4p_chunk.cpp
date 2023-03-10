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

#include "compute_dipole_tip4p_chunk.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "pair.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathSpecial;

enum { MASSCENTER, GEOMCENTER };

/* ---------------------------------------------------------------------- */

ComputeDipoleTIP4PChunk::ComputeDipoleTIP4PChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idchunk(nullptr), massproc(nullptr), masstotal(nullptr), chrgproc(nullptr),
  chrgtotal(nullptr), com(nullptr),
  comall(nullptr), dipole(nullptr), dipoleall(nullptr)
{
  if ((narg != 4) && (narg != 5))
    error->all(FLERR,"Illegal compute dipole/tip4p/chunk command");

  array_flag = 1;
  size_array_cols = 4;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  // ID of compute chunk/atom

  idchunk = utils::strdup(arg[3]);

  usecenter = MASSCENTER;

  if (narg == 5) {
    if (strncmp(arg[4],"geom",4) == 0) usecenter = GEOMCENTER;
    else if (strcmp(arg[4],"mass") == 0) usecenter = MASSCENTER;
    else error->all(FLERR,"Illegal compute dipole/tip4p/chunk command");
  }

  ComputeDipoleTIP4PChunk::init();

  // chunk-based data

  nchunk = 1;
  maxchunk = 0;
  allocate();
}

/* ---------------------------------------------------------------------- */

ComputeDipoleTIP4PChunk::~ComputeDipoleTIP4PChunk()
{
  delete[] idchunk;
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

void ComputeDipoleTIP4PChunk::init()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute < 0)
    error->all(FLERR,"Chunk/atom compute does not exist for compute dipole/tip4p/chunk");
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->compute[icompute]);
  if (strcmp(cchunk->style,"chunk/atom") != 0)
    error->all(FLERR,"Compute dipole/tip4p/chunk does not use chunk/atom compute");

  if (!force->pair) error->all(FLERR, "Pair style must be defined for compute dipole/tip4p/chunk");

  int itmp;
  double *p_qdist = (double *) force->pair->extract("qdist", itmp);
  int *p_typeO = (int *) force->pair->extract("typeO", itmp);
  int *p_typeH = (int *) force->pair->extract("typeH", itmp);
  int *p_typeA = (int *) force->pair->extract("typeA", itmp);
  int *p_typeB = (int *) force->pair->extract("typeB", itmp);
  if (!p_qdist || !p_typeO || !p_typeH || !p_typeA || !p_typeB)
    error->all(FLERR, "Pair style is incompatible with compute dipole/tip4p/chunk");
  typeO = *p_typeO;
  typeH = *p_typeH;
  int typeA = *p_typeA;
  int typeB = *p_typeB;

  if (!force->angle || !force->bond || !force->angle->setflag || !force->bond->setflag)
    error->all(FLERR, "Bond and angle potentials must be defined for compute dipole/tip4p/chunk");
  if ((typeA < 1) || (typeA > atom->nangletypes) || (force->angle->setflag[typeA] == 0))
    error->all(FLERR, "Bad TIP4P angle type for compute dipole/tip4p/chunk");
  if ((typeB < 1) || (typeB > atom->nbondtypes) || (force->bond->setflag[typeB] == 0))
    error->all(FLERR, "Bad TIP4P bond type for compute dipole/tip4p/chunk");
  double theta = force->angle->equilibrium_angle(typeA);
  double blen = force->bond->equilibrium_distance(typeB);
  alpha = *p_qdist / (cos(0.5 * theta) * blen);
}

/* ---------------------------------------------------------------------- */

void ComputeDipoleTIP4PChunk::compute_array()
{
  int i,index;
  double massone;
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
  double xM[3];
  double *xi;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i]-1;
      if (index < 0) continue;
      if (usecenter == MASSCENTER) {
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
      } else massone = 1.0;     // usecenter == GEOMCENTER
      domain->unmap(x[i],image[i],unwrap);
      massproc[index] += massone;
      if (atom->q_flag) chrgproc[index] += q[i];
      com[index][0] += unwrap[0] * massone;
      com[index][1] += unwrap[1] * massone;
      com[index][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(massproc,masstotal,nchunk,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(chrgproc,chrgtotal,nchunk,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&com[0][0],&comall[0][0],3*nchunk,MPI_DOUBLE,MPI_SUM,world);

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
      index = ichunk[i]-1;
      if (index < 0) continue;

      if (type[i] == typeO) {
        find_M(i,xM);
        xi = xM;
      } else {
        xi = x[i];
      }
      domain->unmap(xi,image[i],unwrap);
      if (atom->q_flag) {
        dipole[index][0] += q[i]*unwrap[0];
        dipole[index][1] += q[i]*unwrap[1];
        dipole[index][2] += q[i]*unwrap[2];
      }
      if (atom->mu_flag) {
        dipole[index][0] += mu[i][0];
        dipole[index][1] += mu[i][1];
        dipole[index][2] += mu[i][2];
      }
    }
  }

  MPI_Allreduce(&dipole[0][0],&dipoleall[0][0],4*nchunk,MPI_DOUBLE,MPI_SUM,world);

  for (i = 0; i < nchunk; i++) {
    // correct for position dependence with charged chunks
    dipoleall[i][0] -= chrgtotal[i]*comall[i][0];
    dipoleall[i][1] -= chrgtotal[i]*comall[i][1];
    dipoleall[i][2] -= chrgtotal[i]*comall[i][2];
    // compute total dipole moment
    dipoleall[i][3] = sqrt(square(dipoleall[i][0])
                           + square(dipoleall[i][1])
                           + square(dipoleall[i][2]));
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

void ComputeDipoleTIP4PChunk::lock_enable()
{
  cchunk->lockcount++;
}

/* ----------------------------------------------------------------------
   decrement lock counter in compute chunk/atom, it if still exists
------------------------------------------------------------------------- */

void ComputeDipoleTIP4PChunk::lock_disable()
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

int ComputeDipoleTIP4PChunk::lock_length()
{
  nchunk = cchunk->setup_chunks();
  return nchunk;
}

/* ----------------------------------------------------------------------
   set the lock from startstep to stopstep
------------------------------------------------------------------------- */

void ComputeDipoleTIP4PChunk::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  cchunk->lock(fixptr,startstep,stopstep);
}

/* ----------------------------------------------------------------------
   unset the lock
------------------------------------------------------------------------- */

void ComputeDipoleTIP4PChunk::unlock(Fix *fixptr)
{
  cchunk->unlock(fixptr);
}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeDipoleTIP4PChunk::allocate()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(chrgproc);
  memory->destroy(chrgtotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(dipole);
  memory->destroy(dipoleall);
  maxchunk = nchunk;
  memory->create(massproc,maxchunk,"dipole/tip4p/chunk:massproc");
  memory->create(masstotal,maxchunk,"dipole/tip4p/chunk:masstotal");
  memory->create(chrgproc,maxchunk,"dipole/tip4p/chunk:chrgproc");
  memory->create(chrgtotal,maxchunk,"dipole/tip4p/chunk:chrgtotal");
  memory->create(com,maxchunk,3,"dipole/tip4p/chunk:com");
  memory->create(comall,maxchunk,3,"dipole/tip4p/chunk:comall");
  memory->create(dipole,maxchunk,4,"dipole/tip4p/chunk:dipole");
  memory->create(dipoleall,maxchunk,4,"dipole/tip4p/chunk:dipoleall");
  array = dipoleall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeDipoleTIP4PChunk::memory_usage()
{
  double bytes = (bigint) maxchunk * 2 * sizeof(double);
  bytes += (double)maxchunk * 2*3 * sizeof(double);
  bytes += (double)maxchunk * 2*4 * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void ComputeDipoleTIP4PChunk::find_M(int i, double *xM)
{
  double **x = atom->x;

  int iH1 = atom->map(atom->tag[i] + 1);
  int iH2 = atom->map(atom->tag[i] + 2);

  if ((iH1 == -1) || (iH2 == -1)) error->one(FLERR,"TIP4P hydrogen is missing");
  if ((atom->type[iH1] != typeH) || (atom->type[iH2] != typeH))
    error->one(FLERR,"TIP4P hydrogen has incorrect atom type");

  // set iH1,iH2 to index of closest image to O

  iH1 = domain->closest_image(i,iH1);
  iH2 = domain->closest_image(i,iH2);

  double delx1 = x[iH1][0] - x[i][0];
  double dely1 = x[iH1][1] - x[i][1];
  double delz1 = x[iH1][2] - x[i][2];

  double delx2 = x[iH2][0] - x[i][0];
  double dely2 = x[iH2][1] - x[i][1];
  double delz2 = x[iH2][2] - x[i][2];

  xM[0] = x[i][0] + alpha * 0.5 * (delx1 + delx2);
  xM[1] = x[i][1] + alpha * 0.5 * (dely1 + dely2);
  xM[2] = x[i][2] + alpha * 0.5 * (delz1 + delz2);
}

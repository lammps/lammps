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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_spring_chunk.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "domain.h"
#include "modify.h"
#include "compute_chunk_atom.h"
#include "compute_com_chunk.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1.0e-10

/* ---------------------------------------------------------------------- */

FixSpringChunk::FixSpringChunk(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  idchunk(NULL), idcom(NULL), com0(NULL), fcom(NULL)
{
  if (narg != 6) error->all(FLERR,"Illegal fix spring/chunk command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  k_spring = force->numeric(FLERR,arg[3]);

  int n = strlen(arg[4]) + 1;
  idchunk = new char[n];
  strcpy(idchunk,arg[4]);

  n = strlen(arg[5]) + 1;
  idcom = new char[n];
  strcpy(idcom,arg[5]);

  esprings = 0.0;
  nchunk = 0;
}

/* ---------------------------------------------------------------------- */

FixSpringChunk::~FixSpringChunk()
{
  memory->destroy(com0);
  memory->destroy(fcom);

  // decrement lock counter in compute chunk/atom, it if still exists

  int icompute = modify->find_compute(idchunk);
  if (icompute >= 0) {
    cchunk = (ComputeChunkAtom *) modify->compute[icompute];
    cchunk->unlock(this);
    cchunk->lockcount--;
  }

  delete [] idchunk;
  delete [] idcom;
}

/* ---------------------------------------------------------------------- */

int FixSpringChunk::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpringChunk::init()
{
  // current indices for idchunk and idcom

  int icompute = modify->find_compute(idchunk);
  if (icompute < 0)
    error->all(FLERR,"Chunk/atom compute does not exist for fix spring/chunk");
  cchunk = (ComputeChunkAtom *) modify->compute[icompute];
  if (strcmp(cchunk->style,"chunk/atom") != 0)
    error->all(FLERR,"Fix spring/chunk does not use chunk/atom compute");

  icompute = modify->find_compute(idcom);
  if (icompute < 0)
    error->all(FLERR,"Com/chunk compute does not exist for fix spring/chunk");
  ccom = (ComputeCOMChunk *) modify->compute[icompute];
  if (strcmp(ccom->style,"com/chunk") != 0)
    error->all(FLERR,"Fix spring/chunk does not use com/chunk compute");

  // check that idchunk is consistent with ccom->idchunk

  if (strcmp(idchunk,ccom->idchunk) != 0)
    error->all(FLERR,"Fix spring chunk chunkID not same as comID chunkID");

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringChunk::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
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

  if (com0 == NULL) cchunk->lock(this,update->ntimestep,-1);

  // calculate current centers of mass for each chunk
  // extract pointers from idchunk and idcom

  ccom->compute_array();

  nchunk = cchunk->nchunk;
  int *ichunk = cchunk->ichunk;
  double *masstotal = ccom->masstotal;
  double **com = ccom->array;

  // check if first time cchunk was queried via ccom
  // if so, allocate com0,fcom and store initial COM

  if (com0 == NULL) {
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
   energy of stretched spring
------------------------------------------------------------------------- */

double FixSpringChunk::compute_scalar()
{
  return esprings;
}

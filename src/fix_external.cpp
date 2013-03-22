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

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_external.h"
#include "atom.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{PF_CALLBACK,PF_ARRAY};

/* ---------------------------------------------------------------------- */

FixExternal::FixExternal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix external command");

  if (strcmp(arg[3],"pf/callback") == 0) {
    if (narg != 6) error->all(FLERR,"Illegal fix external command");
    mode = PF_CALLBACK;
    ncall = atoi(arg[4]);
    napply = atoi(arg[5]);
    if (ncall <= 0 || napply <= 0) 
      error->all(FLERR,"Illegal fix external command");
  } else if (strcmp(arg[3],"pf/array") == 0) {
    if (narg != 5) error->all(FLERR,"Illegal fix external command");
    mode = PF_ARRAY;
    napply = atoi(arg[4]);
    if (napply <= 0) error->all(FLERR,"Illegal fix external command");
  } else error->all(FLERR,"Illegal fix external command");

  callback = NULL;

  // perform initial allocation of atom-based array
  // register with Atom class

  fexternal = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
}

/* ---------------------------------------------------------------------- */

FixExternal::~FixExternal()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  memory->destroy(fexternal);
}

/* ---------------------------------------------------------------------- */

int FixExternal::setmask()
{
  int mask = 0;
  if (mode == PF_CALLBACK || mode == PF_ARRAY) {
    mask |= POST_FORCE;
    mask |= MIN_POST_FORCE;
  }
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExternal::init()
{
  if (mode == PF_CALLBACK && callback == NULL)
    error->all(FLERR,"Fix external callback function not set");
}

/* ---------------------------------------------------------------------- */

void FixExternal::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixExternal::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixExternal::post_force(int vflag)
{
  bigint ntimestep = update->ntimestep;

  // invoke the callback in driver program
  // it will fill fexternal with forces

  if (mode == PF_CALLBACK && ntimestep % ncall == 0)
    (this->callback)(ptr_caller,update->ntimestep,
                     atom->nlocal,atom->tag,atom->x,fexternal);

  // add forces from fexternal to atoms in group

  if (ntimestep % napply == 0) {
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        f[i][0] += fexternal[i][0];
        f[i][1] += fexternal[i][1];
        f[i][2] += fexternal[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixExternal::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixExternal::memory_usage()
{
  double bytes = 3*atom->nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixExternal::grow_arrays(int nmax)
{
  memory->grow(fexternal,nmax,3,"external:fexternal");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixExternal::copy_arrays(int i, int j, int delflag)
{
  fexternal[j][0] = fexternal[i][0];
  fexternal[j][1] = fexternal[i][1];
  fexternal[j][2] = fexternal[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixExternal::pack_exchange(int i, double *buf)
{
  buf[0] = fexternal[i][0];
  buf[1] = fexternal[i][1];
  buf[2] = fexternal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixExternal::unpack_exchange(int nlocal, double *buf)
{
  fexternal[nlocal][0] = buf[0];
  fexternal[nlocal][1] = buf[1];
  fexternal[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   external caller sets a callback function to invoke in post_force()
------------------------------------------------------------------------- */

void FixExternal::set_callback(FnPtr caller_callback, void *caller_ptr)
{
  callback = caller_callback;
  ptr_caller = caller_ptr;
}

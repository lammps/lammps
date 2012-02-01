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
#include "fix_external.h"
#include "atom.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixExternal::FixExternal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix external command");

  callback = NULL;

  nmax = 0;
  fexternal = NULL;
}

/* ---------------------------------------------------------------------- */

FixExternal::~FixExternal()
{
  memory->destroy(fexternal);
}

/* ---------------------------------------------------------------------- */

int FixExternal::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExternal::init()
{
  if (callback == NULL) error->all(FLERR,"Fix external callback function not set");
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
  if (atom->nlocal > nmax) {
    memory->destroy(fexternal);
    nmax = atom->nmax;
    memory->create(fexternal,nmax,3,"external:fexternal");
  }

  // invoke the callback in driver program
  // it will fill fexternal with forces

  (this->callback)(ptr_caller,update->ntimestep,
		   atom->nlocal,atom->tag,atom->x,fexternal);

  // add forces from fexternal to atoms in group

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

/* ---------------------------------------------------------------------- */

void FixExternal::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   external caller sets a callback function to invoke in post_force()
------------------------------------------------------------------------- */

void FixExternal::set_callback(FnPtr caller_callback, void *caller_ptr)
{
  callback = caller_callback;
  ptr_caller = caller_ptr;
}

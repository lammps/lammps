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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_wiggle.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWiggle::FixWiggle(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all("Illegal fix wiggle command");

  time_depend = 1;

  // parse command-line args

  if (strcmp(arg[3],"x") == 0) axis = 0;
  else if (strcmp(arg[3],"y") == 0) axis = 1;
  else if (strcmp(arg[3],"z") == 0) axis = 2;
  else error->all("Illegal fix wiggle command");
  
  amplitude = atof(arg[4]);
  period = atof(arg[5]);

  // perform initial allocation of atom-based array
  // register with Atom class
  
  original = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // oscillation parameters

  double PI = 4.0 * atan(1.0);
  omega = 2.0*PI / period;
  time_origin = update->ntimestep;

  // original = initial atom coord in wiggled direction
  
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) original[i] = x[i][axis];
    else original[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixWiggle::~FixWiggle()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored array

  memory->sfree(original);
}

/* ---------------------------------------------------------------------- */

int FixWiggle::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWiggle::init()
{
  dt = update->dt;

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixWiggle::post_force(int vflag)
{
  double arg = omega * (update->ntimestep - time_origin) * dt;
  double cosine = cos(arg);
  double sine = sin(arg);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][axis] = original[i] + amplitude - amplitude*cosine;
      v[i][axis] = amplitude*omega*sine;
      f[i][axis] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWiggle::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array 
------------------------------------------------------------------------- */

double FixWiggle::memory_usage()
{
  double bytes = atom->nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based array 
------------------------------------------------------------------------- */

void FixWiggle::grow_arrays(int nmax)
{
  original = (double *)
    memory->srealloc(original,nmax*sizeof(double),"wiggle:original");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array 
------------------------------------------------------------------------- */

void FixWiggle::copy_arrays(int i, int j)
{
  original[j] = original[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc 
------------------------------------------------------------------------- */

int FixWiggle::pack_exchange(int i, double *buf)
{
  buf[0] = original[i];
  return 1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc 
------------------------------------------------------------------------- */

int FixWiggle::unpack_exchange(int nlocal, double *buf)
{
  original[nlocal] = buf[0];
  return 1;
}

/* ---------------------------------------------------------------------- */

void FixWiggle::reset_dt()
{
  dt = update->dt;
}

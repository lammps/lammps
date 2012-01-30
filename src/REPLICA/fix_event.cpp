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

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL), Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "fix_event.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "neighbor.h"
#include "comm.h"
#include "universe.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEvent::FixEvent(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix event command");

  restart_global = 1;

  // perform initial allocation of atom-based array
  // register with Atom class

  xevent = NULL;
  xold = NULL;
  vold = NULL;
  imageold = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
}

/* ---------------------------------------------------------------------- */

FixEvent::~FixEvent()
{
  // unregister callbacks to this fix from Atom class
 
  atom->delete_callback(id,0);

  // delete locally stored array

  memory->destroy(xevent);
  memory->destroy(xold);
  memory->destroy(vold);
  memory->destroy(imageold);
}

/* ---------------------------------------------------------------------- */

int FixEvent::setmask()
{
  return 0;
}

/* ----------------------------------------------------------------------
   save current atom coords as an event
   called when an event occurs
------------------------------------------------------------------------- */

void FixEvent::store_event()
{
  double **x = atom->x;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) 
    domain->unmap(x[i],image[i],xevent[i]);

}

/* ----------------------------------------------------------------------
   restore atom coords to quenched initial state
   called prior to NEB calculation
------------------------------------------------------------------------- */

void FixEvent::restore_event()
{
  double **x = atom->x;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    x[i][0] = xevent[i][0];
    x[i][1] = xevent[i][1];
    x[i][2] = xevent[i][2];

    // Since xevent is unwrapped coordinate, need to
    // adjust image flags when remapping

    image[i] = (512 << 20) | (512 << 10) | 512;
    domain->remap(x[i],image[i]);
    //    domain->remap(x[i]);
  }

}

/* ----------------------------------------------------------------------
   store state of all atoms
   called before quench and subsequent check for event
   so can later restore pre-quench state if no event occurs
------------------------------------------------------------------------- */

void FixEvent::store_state()
{
  double **x = atom->x;
  double **v = atom->v;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    xold[i][0] = x[i][0];
    xold[i][1] = x[i][1];
    xold[i][2] = x[i][2];
    vold[i][0] = v[i][0];
    vold[i][1] = v[i][1];
    vold[i][2] = v[i][2];
    imageold[i] = image[i];
  }
}

/* ----------------------------------------------------------------------
   restore state of all atoms to pre-quench state
   called after no event detected so can continue
------------------------------------------------------------------------- */

void FixEvent::restore_state()
{
  double **x = atom->x;
  double **v = atom->v;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    x[i][0] = xold[i][0];
    x[i][1] = xold[i][1];
    x[i][2] = xold[i][2];
    v[i][0] = vold[i][0];
    v[i][1] = vold[i][1];
    v[i][2] = vold[i][2];
    image[i] = imageold[i];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixEvent::memory_usage()
{
  double bytes = 6*atom->nmax * sizeof(double);
  bytes += atom->nmax*sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixEvent::grow_arrays(int nmax)
{
  memory->grow(xevent,nmax,3,"event:xevent");
  memory->grow(xold,nmax,3,"event:xold");
  memory->grow(vold,nmax,3,"event:vold");
  memory->grow(imageold,nmax,"event:imageold");

  // allow compute event to access stored event coords

  array_atom = xevent;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixEvent::copy_arrays(int i, int j)
{
  xevent[j][0] = xevent[i][0];
  xevent[j][1] = xevent[i][1];
  xevent[j][2] = xevent[i][2];
  xold[j][0] = xold[i][0];
  xold[j][1] = xold[i][1];
  xold[j][2] = xold[i][2];
  vold[j][0] = vold[i][0];
  vold[j][1] = vold[i][1];
  vold[j][2] = vold[i][2];
  imageold[j] = imageold[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixEvent::pack_exchange(int i, double *buf)
{
  buf[0] = xevent[i][0];
  buf[1] = xevent[i][1];
  buf[2] = xevent[i][2];
  buf[3] = xold[i][0];
  buf[4] = xold[i][1];
  buf[5] = xold[i][2];
  buf[6] = vold[i][0];
  buf[7] = vold[i][1];
  buf[8] = vold[i][2];
  buf[9] = imageold[i];

  return 10;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixEvent::unpack_exchange(int nlocal, double *buf)
{
  xevent[nlocal][0] = buf[0];
  xevent[nlocal][1] = buf[1];
  xevent[nlocal][2] = buf[2];
  xold[nlocal][0] = buf[3];
  xold[nlocal][1] = buf[4];
  xold[nlocal][2] = buf[5];
  vold[nlocal][0] = buf[6];
  vold[nlocal][1] = buf[7];
  vold[nlocal][2] = buf[8];
  imageold[nlocal] = static_cast<int>(buf[9]);

  return 10;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write 
------------------------------------------------------------------------- */

void FixEvent::write_restart(FILE *fp)
{
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix 
------------------------------------------------------------------------- */

void FixEvent::restart(char *buf)
{
}

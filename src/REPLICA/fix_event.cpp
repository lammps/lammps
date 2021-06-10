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

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL), Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "fix_event.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEvent::FixEvent(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), xevent(nullptr), xold(nullptr), vold(nullptr),
  imageold(nullptr), xorig(nullptr), vorig(nullptr), imageorig(nullptr)
{
  if (narg != 3) error->all(FLERR,"Illegal fix event command");

  restart_global = 1;

  // perform initial allocation of atom-based array
  // register with Atom class

  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
}

/* ---------------------------------------------------------------------- */

FixEvent::~FixEvent()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);

  // delete locally stored array

  memory->destroy(xevent);
  memory->destroy(xold);
  memory->destroy(vold);
  memory->destroy(imageold);
  memory->destroy(xorig);
  memory->destroy(vorig);
  memory->destroy(imageorig);
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
  imageint *image = atom->image;
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
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    x[i][0] = xevent[i][0];
    x[i][1] = xevent[i][1];
    x[i][2] = xevent[i][2];

    // since xevent is unwrapped coordinate,
    // need to adjust image flag from default when remapping
    // same as in read_data -> Atom::data_atoms()

    image[i] = ((imageint) IMGMAX << IMG2BITS) |
      ((imageint) IMGMAX << IMGBITS) | IMGMAX;
    domain->remap(x[i],image[i]);
  }

}

/* ----------------------------------------------------------------------
   store state of all atoms
   called before quench and subsequent check for event
   so can later restore pre-quench state if no event occurs
------------------------------------------------------------------------- */

void FixEvent::store_state_quench()
{
  double **x = atom->x;
  double **v = atom->v;
  imageint *image = atom->image;
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

void FixEvent::restore_state_quench()
{
  double **x = atom->x;
  double **v = atom->v;
  imageint *image = atom->image;
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
   store original state of all atoms
------------------------------------------------------------------------- */

void FixEvent::store_state_dephase()
{
  double **x = atom->x;
  double **v = atom->v;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    xorig[i][0] = x[i][0];
    xorig[i][1] = x[i][1];
    xorig[i][2] = x[i][2];
    vorig[i][0] = v[i][0];
    vorig[i][1] = v[i][1];
    vorig[i][2] = v[i][2];
    imageorig[i] = image[i];
  }
}

/* ----------------------------------------------------------------------
   restore state of all atoms to original state
------------------------------------------------------------------------- */

void FixEvent::restore_state_dephase()
{
  double **x = atom->x;
  double **v = atom->v;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    x[i][0] = xorig[i][0];
    x[i][1] = xorig[i][1];
    x[i][2] = xorig[i][2];
    v[i][0] = vorig[i][0];
    v[i][1] = vorig[i][1];
    v[i][2] = vorig[i][2];
    image[i] = imageorig[i];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixEvent::memory_usage()
{
  double bytes = 12*atom->nmax * sizeof(double);
  bytes += (double)atom->nmax*sizeof(int);
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
  memory->grow(xorig,nmax,3,"event:xorig");
  memory->grow(vorig,nmax,3,"event:vorig");
  memory->grow(imageorig,nmax,"event:imageorig");

  // allow compute event to access stored event coords

  array_atom = xevent;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixEvent::copy_arrays(int i, int j, int /*delflag*/)
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
  xorig[j][0] = xorig[i][0];
  xorig[j][1] = xorig[i][1];
  xorig[j][2] = xorig[i][2];
  vorig[j][0] = vorig[i][0];
  vorig[j][1] = vorig[i][1];
  vorig[j][2] = vorig[i][2];
  imageorig[j] = imageorig[i];
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
  buf[10] = xorig[i][0];
  buf[11] = xorig[i][1];
  buf[12] = xorig[i][2];
  buf[13] = vorig[i][0];
  buf[14] = vorig[i][1];
  buf[15] = vorig[i][2];
  buf[16] = imageorig[i];

  return 17;
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
  imageold[nlocal] = static_cast<imageint>(buf[9]);
  xorig[nlocal][0] = buf[10];
  xorig[nlocal][1] = buf[11];
  xorig[nlocal][2] = buf[12];
  vorig[nlocal][0] = buf[13];
  vorig[nlocal][1] = buf[14];
  vorig[nlocal][2] = buf[15];
  imageorig[nlocal] = static_cast<imageint>(buf[16]);

  return 17;
}

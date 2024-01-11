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

#include "fix_store_force.h"

#include "atom.h"
#include "error.h"
#include "memory.h"
#include "respa.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixStoreForce::FixStoreForce(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), foriginal(nullptr)
{
  if (narg < 3) error->all(FLERR, "Illegal fix store/force command");

  peratom_flag = 1;
  size_peratom_cols = 3;
  peratom_freq = 1;

  nmax = atom->nmax;
  memory->create(foriginal, nmax, 3, "store/force:foriginal");
  array_atom = foriginal;

  // zero the array since dump may access it on timestep 0
  // zero the array since a variable may access it before first run

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) foriginal[i][0] = foriginal[i][1] = foriginal[i][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixStoreForce::~FixStoreForce()
{
  memory->destroy(foriginal);
}

/* ---------------------------------------------------------------------- */

int FixStoreForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStoreForce::init()
{
  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixStoreForce::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa - 1);
    post_force_respa(vflag, nlevels_respa - 1, 0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa - 1);
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixStoreForce::post_force(int /*vflag*/)
{
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(foriginal);
    memory->create(foriginal, nmax, 3, "store/force:foriginal");
    array_atom = foriginal;
  }

  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      foriginal[i][0] = f[i][0];
      foriginal[i][1] = f[i][1];
      foriginal[i][2] = f[i][2];
    } else
      foriginal[i][0] = foriginal[i][1] = foriginal[i][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixStoreForce::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixStoreForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixStoreForce::memory_usage()
{
  double bytes = (double) atom->nmax * 3 * sizeof(double);
  return bytes;
}

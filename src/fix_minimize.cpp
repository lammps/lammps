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

#include "fix_minimize.h"
#include "atom.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixMinimize::FixMinimize(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // perform initial allocation of atom-based arrays
  // register with Atom class

  gradient = NULL;
  searchdir = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
}

/* ---------------------------------------------------------------------- */

FixMinimize::~FixMinimize()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored arrays

  memory->destroy_2d_double_array(gradient);
  memory->destroy_2d_double_array(searchdir);
}

/* ---------------------------------------------------------------------- */

int FixMinimize::setmask()
{
  return 0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double FixMinimize::memory_usage()
{
  double bytes = 2 * atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays 
------------------------------------------------------------------------- */

void FixMinimize::grow_arrays(int nmax)
{
  gradient =
    memory->grow_2d_double_array(gradient,nmax,3,"fix_minimize:gradient");
  searchdir =
    memory->grow_2d_double_array(searchdir,nmax,3,"fix_minimize:searchdir");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays 
------------------------------------------------------------------------- */

void FixMinimize::copy_arrays(int i, int j)
{
  gradient[j][0] = gradient[i][0];
  gradient[j][1] = gradient[i][1];
  gradient[j][2] = gradient[i][2];
  searchdir[j][0] = searchdir[i][0];
  searchdir[j][1] = searchdir[i][1];
  searchdir[j][2] = searchdir[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc 
------------------------------------------------------------------------- */

int FixMinimize::pack_exchange(int i, double *buf)
{
  buf[0] = gradient[i][0]; 
  buf[1] = gradient[i][1]; 
  buf[2] = gradient[i][2]; 
  buf[3] = searchdir[i][0]; 
  buf[4] = searchdir[i][1]; 
  buf[5] = searchdir[i][2]; 
  return 6;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc 
------------------------------------------------------------------------- */

int FixMinimize::unpack_exchange(int nlocal, double *buf)
{
  gradient[nlocal][0] = buf[0];
  gradient[nlocal][1] = buf[1];
  gradient[nlocal][2] = buf[2];
  searchdir[nlocal][0] = buf[3];
  searchdir[nlocal][1] = buf[4];
  searchdir[nlocal][2] = buf[5];
  return 6;
}

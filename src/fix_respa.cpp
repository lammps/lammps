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

#include "stdlib.h"
#include "fix_respa.h"
#include "atom.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixRespa::FixRespa(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // nlevels = # of rRESPA levels

  nlevels = atoi(arg[3]);

  // perform initial allocation of atom-based arrays
  // register with Atom class

  f_level = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
}

/* ---------------------------------------------------------------------- */

FixRespa::~FixRespa()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored arrays

  memory->destroy_3d_double_array(f_level);
}

/* ---------------------------------------------------------------------- */

int FixRespa::setmask()
{
  return 0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double FixRespa::memory_usage()
{
  double bytes = atom->nmax*nlevels*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays 
------------------------------------------------------------------------- */

void FixRespa::grow_arrays(int nmax)
{
  f_level = 
    memory->grow_3d_double_array(f_level,nmax,nlevels,3,"fix_respa:f_level");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays 
------------------------------------------------------------------------- */

void FixRespa::copy_arrays(int i, int j)
{
  for (int k = 0; k < nlevels; k++) {
    f_level[j][k][0] = f_level[i][k][0];
    f_level[j][k][1] = f_level[i][k][1];
    f_level[j][k][2] = f_level[i][k][2];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc 
------------------------------------------------------------------------- */

int FixRespa::pack_exchange(int i, double *buf)
{
  int m = 0;
  for (int k = 0; k < nlevels; k++) {
    buf[m++] = f_level[i][k][0];
    buf[m++] = f_level[i][k][1];
    buf[m++] = f_level[i][k][2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc 
------------------------------------------------------------------------- */

int FixRespa::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  for (int k = 0; k < nlevels; k++) {
    f_level[nlocal][k][0] = buf[m++];
    f_level[nlocal][k][1] = buf[m++];
    f_level[nlocal][k][2] = buf[m++];
  }
  return m;
}

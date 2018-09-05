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

#include <cstdlib>
#include "fix_read_restart.h"
#include "atom.h"
#include "memory.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReadRestart::FixReadRestart(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  count(NULL), extra(NULL)
{
  nextra = force->inumeric(FLERR,arg[3]);
  int nfix = force->inumeric(FLERR,arg[4]);

  // perform initial allocation of atom-based array
  // register with Atom class

  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // extra = copy of atom->extra

  double **atom_extra = atom->extra;
  int nlocal = atom->nlocal;
  int i,j,m;

  for (i = 0; i < nlocal; i++) {
    m = 0;
    for (j = 0; j < nfix; j++) m += static_cast<int> (atom_extra[i][m]);
    count[i] = m;
    for (j = 0; j < m; j++) extra[i][j] = atom_extra[i][j];
  }
}

/* ---------------------------------------------------------------------- */

FixReadRestart::~FixReadRestart()
{
  // unregister callback to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored arrays

  memory->destroy(count);
  memory->destroy(extra);
}

/* ---------------------------------------------------------------------- */

int FixReadRestart::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixReadRestart::memory_usage()
{
  double bytes = atom->nmax*nextra * sizeof(double);
  bytes += atom->nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixReadRestart::grow_arrays(int nmax)
{
  memory->grow(count,nmax,"read_restart:count");
  memory->grow(extra,nmax,nextra,"read_restart:extra");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixReadRestart::copy_arrays(int i, int j, int /*delflag*/)
{
  count[j] = count[i];
  for (int m = 0; m < count[i]; m++) extra[j][m] = extra[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixReadRestart::pack_exchange(int i, double *buf)
{
  buf[0] = count[i];
  for (int m = 0; m < count[i]; m++) buf[m+1] = extra[i][m];
  return count[i]+1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixReadRestart::unpack_exchange(int nlocal, double *buf)
{
  count[nlocal] = static_cast<int> (buf[0]);
  for (int m = 0; m < count[nlocal]; m++) extra[nlocal][m] = buf[m+1];
  return count[nlocal]+1;
}

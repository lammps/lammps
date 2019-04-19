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
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Please cite the related publication:
   H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
   "Parallel Reactive Molecular Dynamics: Numerical Methods and
   Algorithmic Techniques", Parallel Computing, in press.
------------------------------------------------------------------------- */

#include "fix_reaxc.h"
#include "atom.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAX_REAX_BONDS      30
#define MIN_REAX_BONDS      15
#define MIN_REAX_HBONDS     25

/* ---------------------------------------------------------------------- */

FixReaxC::FixReaxC(LAMMPS *lmp,int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // perform initial allocation of atom-based arrays
  // register with atom class

  oldnmax = 0;
  num_bonds = NULL;
  num_hbonds = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // initialize arrays to MIN so atom migration is OK the 1st time
  // it is done in grow_arrays() now

  // set comm sizes needed by this fix

  comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixReaxC::~FixReaxC()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);

  // delete locally stored arrays

  memory->destroy(num_bonds);
  memory->destroy(num_hbonds);
}

/* ---------------------------------------------------------------------- */

int FixReaxC::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixReaxC::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * 2 * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixReaxC::grow_arrays(int nmax)
{
  memory->grow(num_bonds,nmax,"reaxc:num_bonds");
  memory->grow(num_hbonds,nmax,"reaxc:num_hbonds");
  for (int i = oldnmax; i < nmax; i++) {
    num_hbonds[i] = MIN_REAX_HBONDS;
    num_bonds[i] = MIN_REAX_BONDS;
  }
  oldnmax = nmax;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixReaxC::copy_arrays(int i, int j, int /*delflag*/)
{
  num_bonds[j] = num_bonds[i];
  num_hbonds[j] = num_hbonds[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixReaxC::pack_exchange(int i, double *buf)
{
  buf[0] = num_bonds[i];
  buf[1] = num_hbonds[i];
  return 2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixReaxC::unpack_exchange(int nlocal, double *buf)
{
  num_bonds[nlocal] = static_cast<int> (buf[0]);
  num_hbonds[nlocal] = static_cast<int> (buf[1]);
  return 2;
}

/* ---------------------------------------------------------------------- */

int FixReaxC::pack_forward_comm(int n, int *list, double *buf,
                                int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = num_bonds[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixReaxC::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    num_bonds[i] = static_cast<int> (buf[m++]);
}

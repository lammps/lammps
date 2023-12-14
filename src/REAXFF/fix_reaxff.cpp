// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Please cite the related publication:
   H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
   "Parallel Reactive Molecular Dynamics: Numerical Methods and
   Algorithmic Techniques", Parallel Computing, in press.
------------------------------------------------------------------------- */

#include "fix_reaxff.h"
#include "atom.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAX_REAX_BONDS      30
#define MIN_REAX_BONDS      15
#define MIN_REAX_HBONDS     25

/* ---------------------------------------------------------------------- */

FixReaxFF::FixReaxFF(LAMMPS *lmp,int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // perform initial allocation of atom-based arrays
  // register with atom class

  oldnmax = 0;
  num_bonds = nullptr;
  num_hbonds = nullptr;
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // initialize arrays to MIN so atom migration is OK the 1st time
  // it is done in grow_arrays() now

  // set comm sizes needed by this fix

  comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixReaxFF::~FixReaxFF()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,Atom::GROW);

  // delete locally stored arrays

  memory->destroy(num_bonds);
  memory->destroy(num_hbonds);
}

/* ---------------------------------------------------------------------- */

int FixReaxFF::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixReaxFF::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = (double)nmax * 2 * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixReaxFF::grow_arrays(int nmax)
{
  memory->grow(num_bonds,nmax,"reaxff:num_bonds");
  memory->grow(num_hbonds,nmax,"reaxff:num_hbonds");
  for (int i = oldnmax; i < nmax; i++) {
    num_hbonds[i] = MIN_REAX_HBONDS;
    num_bonds[i] = MIN_REAX_BONDS;
  }
  oldnmax = nmax;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixReaxFF::copy_arrays(int i, int j, int /*delflag*/)
{
  num_bonds[j] = num_bonds[i];
  num_hbonds[j] = num_hbonds[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixReaxFF::pack_exchange(int i, double *buf)
{
  buf[0] = num_bonds[i];
  buf[1] = num_hbonds[i];
  return 2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixReaxFF::unpack_exchange(int nlocal, double *buf)
{
  num_bonds[nlocal] = static_cast<int> (buf[0]);
  num_hbonds[nlocal] = static_cast<int> (buf[1]);
  return 2;
}

/* ---------------------------------------------------------------------- */

int FixReaxFF::pack_forward_comm(int n, int *list, double *buf,
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

void FixReaxFF::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    num_bonds[i] = static_cast<int> (buf[m++]);
}

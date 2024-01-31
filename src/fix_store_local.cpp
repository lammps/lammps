/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_store_local.h"

#include "atom.h"
#include "error.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr int DELTA = 1024;

/* ---------------------------------------------------------------------- */

FixStoreLocal::FixStoreLocal(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nvalues(0), vector(nullptr), array(nullptr)
{
  if (narg != 5) error->all(FLERR, "Illegal fix STORE/LOCAL command");
  local_flag = 1;

  nreset = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nreset <= 0) error->all(FLERR, "Illegal fix STORE/LOCAL command");
  local_freq = nreset;

  nvalues = utils::inumeric(FLERR, arg[4], false, lmp);

  if (nvalues <= 0) error->all(FLERR, "Illegal fix STORE/LOCAL command");
  if (nvalues == 1)
    size_local_cols = 0;
  else
    size_local_cols = nvalues;
  size_local_rows = 0;

  vector = nullptr;
  array = nullptr;
  nmax = 0;
  ncount = 0;
}

/* ---------------------------------------------------------------------- */

FixStoreLocal::~FixStoreLocal()
{
  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixStoreLocal::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStoreLocal::add_data(double *input_data, int i, int j)
{
  int *mask = atom->mask;
  if (!(mask[i] & groupbit)) return;
  if (!(mask[j] & groupbit)) return;

  if (ncount >= nmax) reallocate(ncount);

  // fill vector or array with local values
  if (nvalues == 1) {
    vector[ncount] = input_data[0];
  } else {
    for (int n = 0; n < nvalues; n++) array[ncount][n] = input_data[n];
  }

  ncount += 1;
}

/* ---------------------------------------------------------------------- */

void FixStoreLocal::post_force(int /*vflag*/)
{
  if (update->ntimestep % nreset == 0) {
    size_local_rows = ncount;
    ncount = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixStoreLocal::reallocate(int n)
{
  // grow vector or array
  while (nmax <= n) nmax += DELTA;

  if (nvalues == 1) {
    memory->grow(vector, nmax, "fix_store_local:vector");
    vector_local = vector;
  } else {
    memory->grow(array, nmax, nvalues, "fix_store_local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double FixStoreLocal::memory_usage()
{
  double bytes = (double) nmax * (double) nvalues * sizeof(double);
  return bytes;
}

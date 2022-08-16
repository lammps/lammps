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

#include "fix_store_peratom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixStorePeratom::FixStorePeratom(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), vstore(nullptr), astore(nullptr)
{
  if (narg != 5 && narg != 6)
    error->all(FLERR, "Illegal fix STORE/PERATOM command: number of args");

  // syntax: id group style 0/1 n1 n2 (n3), last arg optional
  //   0/1 flag = not-store or store peratom values in restart file
  //   N2 = 1 and no n3 is vector, N2 > 1 and no n3 is array, N3 = tensor
  //   nvalues = # of peratom values, N = 1 is vector, N > 1 is array

  disable = 0;
  vecflag = arrayflag = tensorflag = 0;

  restart_peratom = utils::inumeric(FLERR, arg[3], false, lmp);
  n2 = utils::inumeric(FLERR, arg[4], false, lmp);
  if (narg == 6)
    n3 = utils::inumeric(FLERR, arg[5], false, lmp);
  else
    n3 = 1;
  if (restart_peratom < 0 || restart_peratom > 1)
    error->all(FLERR, "Illegal fix STORE/PERATOM restart flag: {}", restart_peratom);
  if (n2 <= 0 || n3 <= 0)
    error->all(FLERR, "Illegal fix STORE/PERATOM dimension args: must be >0: {} {}", n2, n3);
  if (n2 == 1 && narg == 5)
    vecflag = 1;
  else if (narg == 5)
    arrayflag = 1;
  else
    tensorflag = 1;
  nvalues = n2 * n3;
  nbytes = nvalues * sizeof(double);

  vstore = nullptr;
  astore = nullptr;
  tstore = nullptr;

  // allocate data structs and register with Atom class

  FixStorePeratom::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  if (restart_peratom) atom->add_callback(Atom::RESTART);

  // zero the storage

  int nlocal = atom->nlocal;
  if (vecflag) {
    for (int i = 0; i < nlocal; i++) vstore[i] = 0.0;
  } else if (arrayflag) {
    for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < n2; j++) astore[i][j] = 0.0;
  } else if (tensorflag) {
    for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < n2; j++)
        for (int k = 0; k < n3; k++) tstore[i][j][k] = 0.0;
  }
  maxexchange = nvalues;
}

/* ---------------------------------------------------------------------- */

FixStorePeratom::~FixStorePeratom()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id, Atom::GROW);
  if (restart_peratom) atom->delete_callback(id, Atom::RESTART);

  memory->destroy(vstore);
  memory->destroy(astore);
  memory->destroy(tstore);
}

/* ---------------------------------------------------------------------- */

int FixStorePeratom::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixStorePeratom::grow_arrays(int nmax)
{
  if (vecflag)
    memory->grow(vstore, nmax, "store:vstore");
  else if (arrayflag)
    memory->grow(astore, nmax, n2, "store:astore");
  else if (tensorflag)
    memory->grow(tstore, nmax, n2, n3, "store:tstore");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixStorePeratom::copy_arrays(int i, int j, int /*delflag*/)
{
  if (disable) return;

  if (vecflag) {
    vstore[j] = vstore[i];
  } else if (arrayflag) {
    for (int m = 0; m < nvalues; m++) astore[j][m] = astore[i][m];
  } else if (tensorflag) {
    memcpy(&tstore[j][0][0], &tstore[i][0][0], nbytes);
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixStorePeratom::pack_exchange(int i, double *buf)
{
  if (disable) return 0;

  if (vecflag) {
    buf[0] = vstore[i];
  } else if (arrayflag) {
    for (int m = 0; m < nvalues; m++) buf[m] = astore[i][m];
  } else if (tensorflag) {
    memcpy(buf, &tstore[i][0][0], nbytes);
  }

  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixStorePeratom::unpack_exchange(int nlocal, double *buf)
{
  if (disable) return 0;

  if (vecflag) {
    vstore[nlocal] = buf[0];
  } else if (arrayflag) {
    for (int m = 0; m < nvalues; m++) astore[nlocal][m] = buf[m];
  } else if (tensorflag) {
    memcpy(&tstore[nlocal][0][0], buf, nbytes);
  }

  return nvalues;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixStorePeratom::pack_restart(int i, double *buf)
{
  if (disable) {
    buf[0] = 0;
    return 1;
  }

  // pack buf[0] this way because other fixes unpack it
  buf[0] = nvalues + 1;

  if (vecflag) {
    buf[1] = vstore[i];
  } else if (arrayflag) {
    for (int m = 0; m < nvalues; m++) buf[m + 1] = astore[i][m];
  } else if (tensorflag) {
    memcpy(&buf[1], &tstore[i][0][0], nbytes);
  }

  return nvalues + 1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixStorePeratom::unpack_restart(int nlocal, int nth)
{
  if (disable) return;

  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  m++;

  if (vecflag) {
    vstore[nlocal] = extra[nlocal][m];
  } else if (arrayflag) {
    for (int i = 0; i < nvalues; i++) astore[nlocal][i] = extra[nlocal][m++];
  } else if (tensorflag) {
    memcpy(&tstore[nlocal][0][0], &extra[nlocal][m], nbytes);
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixStorePeratom::maxsize_restart()
{
  if (disable) return 1;
  return nvalues + 1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixStorePeratom::size_restart(int /*nlocal*/)
{
  if (disable) return 1;
  return nvalues + 1;
}

/* ----------------------------------------------------------------------
   memory usage of global or peratom atom-based array
------------------------------------------------------------------------- */

double FixStorePeratom::memory_usage()
{
  return (double) atom->nmax * n2 * n3 * sizeof(double);
}

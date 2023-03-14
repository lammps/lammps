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

#include "fix_store_atom.h"

#include "atom.h"
#include "error.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

// INTERNAL fix for storing/communicating per-atom values
// syntax: id group style n1 n2 gflag rflag
//   N1 = 1, N2 = 0 is per-atom vector, single value per atom
//   N1 > 1, N2 = 0 is per-atom array, N1 values per atom
//   N1 > 0, N2 > 0 is per-atom tensor, N1xN2 array per atom
//   gflag = 0/1, no/yes communicate per-atom values with ghost atoms
//   rflag = 0/1, no/yes store per-atom values in restart file

/* ---------------------------------------------------------------------- */

FixStoreAtom::FixStoreAtom(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), vstore(nullptr), astore(nullptr)
{
  if (narg != 7) error->all(FLERR, "Illegal fix STORE/ATOM command: number of args");

  disable = 0;

  n1 = utils::inumeric(FLERR, arg[3], false, lmp);
  n2 = utils::inumeric(FLERR, arg[4], false, lmp);
  ghostflag = utils::logical(FLERR, arg[5], false, lmp);
  restartflag = utils::logical(FLERR, arg[6], false, lmp);

  vecflag = arrayflag = tensorflag = 0;
  if (n1 == 1 && n2 == 0)
    vecflag = 1;
  else if (n1 > 1 && n2 == 0)
    arrayflag = 1;
  else if (n1 > 0 && n2 > 0)
    tensorflag = 1;
  else
    error->all(FLERR, "Illegal fix STORE/ATOM dimension args: {} {}", n1, n2);

  if (vecflag || arrayflag)
    nvalues = n1;
  else
    nvalues = n1 * n2;
  nbytes = nvalues * sizeof(double);

  if (ghostflag) comm_border = nvalues;
  maxexchange = nvalues;

  if (restartflag) restart_peratom = 1;

  // allocate data structs and register with Atom class

  vstore = nullptr;
  astore = nullptr;
  tstore = nullptr;

  FixStoreAtom::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  if (restartflag) atom->add_callback(Atom::RESTART);
  if (ghostflag) atom->add_callback(Atom::BORDER);

  // zero the storage

  int nlocal = atom->nlocal;
  if (vecflag) {
    for (int i = 0; i < nlocal; i++) vstore[i] = 0.0;
  } else if (arrayflag) {
    for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < n1; j++) astore[i][j] = 0.0;
  } else if (tensorflag) {
    for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < n1; j++)
        for (int k = 0; k < n2; k++) tstore[i][j][k] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixStoreAtom::~FixStoreAtom()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id, Atom::GROW);
  if (restartflag) atom->delete_callback(id, Atom::RESTART);
  if (ghostflag) atom->delete_callback(id, Atom::BORDER);

  memory->destroy(vstore);
  memory->destroy(astore);
  memory->destroy(tstore);
}

/* ---------------------------------------------------------------------- */

int FixStoreAtom::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixStoreAtom::grow_arrays(int nmax)
{
  if (vecflag)
    memory->grow(vstore, nmax, "store:vstore");
  else if (arrayflag)
    memory->grow(astore, nmax, n1, "store:astore");
  else if (tensorflag)
    memory->grow(tstore, nmax, n1, n2, "store:tstore");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixStoreAtom::copy_arrays(int i, int j, int /*delflag*/)
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
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixStoreAtom::pack_border(int n, int *list, double *buf)
{
  int i, j, k;

  int m = 0;
  if (vecflag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = vstore[j];
    }
  } else if (arrayflag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      for (k = 0; k < nvalues; k++) buf[m++] = astore[j][k];
    }
  } else if (tensorflag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      memcpy(&buf[m], &tstore[j][0][0], nbytes);
      m += nvalues;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixStoreAtom::unpack_border(int n, int first, double *buf)
{
  int i, k, last;

  int m = 0;
  last = first + n;
  if (vecflag) {
    for (i = first; i < last; i++) vstore[i] = buf[m++];
  } else if (arrayflag) {
    for (i = first; i < last; i++)
      for (k = 0; k < nvalues; k++) astore[i][k] = buf[m++];
  } else if (tensorflag) {
    for (i = first; i < last; i++) {
      memcpy(&tstore[i][0][0], &buf[m], nbytes);
      m += nvalues;
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixStoreAtom::pack_exchange(int i, double *buf)
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

int FixStoreAtom::unpack_exchange(int nlocal, double *buf)
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

int FixStoreAtom::pack_restart(int i, double *buf)
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

void FixStoreAtom::unpack_restart(int nlocal, int nth)
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

int FixStoreAtom::maxsize_restart()
{
  if (disable) return 1;
  return nvalues + 1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixStoreAtom::size_restart(int /*nlocal*/)
{
  if (disable) return 1;
  return nvalues + 1;
}

/* ----------------------------------------------------------------------
   memory usage of per-atom atom-based array
------------------------------------------------------------------------- */

double FixStoreAtom::memory_usage()
{
  return (double) atom->nmax * nvalues * sizeof(double);
}

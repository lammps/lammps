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

#include "compute_force_tally.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "pair.h"
#include "update.h"
#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeForceTally::ComputeForceTally(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR, "Illegal compute force/tally command");

  igroup2 = group->find(arg[3]);
  if (igroup2 == -1) error->all(FLERR, "Could not find compute force/tally second group ID");
  groupbit2 = group->bitmask[igroup2];

  scalar_flag = 1;
  vector_flag = 0;
  peratom_flag = 1;
  timeflag = 1;

  comm_reverse = size_peratom_cols = 3;
  extscalar = 1;
  peflag = 1;    // we need Pair::ev_tally() to be run

  did_setup = invoked_peratom = invoked_scalar = -1;
  nmax = -1;
  fatom = nullptr;
  vector = new double[size_peratom_cols];
}

/* ---------------------------------------------------------------------- */

ComputeForceTally::~ComputeForceTally()
{
  if (force && force->pair) force->pair->del_tally_callback(this);
  memory->destroy(fatom);
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeForceTally::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Trying to use compute force/tally without pair style");
  else
    force->pair->add_tally_callback(this);

  if (comm->me == 0) {
    if (force->pair->single_enable == 0 || force->pair->manybody_flag)
      error->warning(FLERR, "Compute force/tally used with incompatible pair style");

    if (force->bond || force->angle || force->dihedral || force->improper || force->kspace)
      error->warning(FLERR, "Compute force/tally only called from pair style");
  }
  did_setup = -1;
}

/* ---------------------------------------------------------------------- */

void ComputeForceTally::pair_setup_callback(int, int)
{
  // run setup only once per time step.
  // we may be called from multiple pair styles

  if (did_setup == update->ntimestep) return;

  const int ntotal = atom->nlocal + atom->nghost;

  // grow per-atom storage, if needed

  if (atom->nmax > nmax) {
    memory->destroy(fatom);
    nmax = atom->nmax;
    memory->create(fatom, nmax, size_peratom_cols, "force/tally:fatom");
    array_atom = fatom;
  }

  // clear storage

  for (int i = 0; i < ntotal; ++i)
    for (int j = 0; j < size_peratom_cols; ++j) fatom[i][j] = 0.0;

  for (int i = 0; i < size_peratom_cols; ++i) vector[i] = ftotal[i] = 0.0;

  did_setup = update->ntimestep;
}

/* ---------------------------------------------------------------------- */
void ComputeForceTally::pair_tally_callback(int i, int j, int nlocal, int newton, double, double,
                                            double fpair, double dx, double dy, double dz)
{
  const int *const mask = atom->mask;

  if (((mask[i] & groupbit) && (mask[j] & groupbit2)) ||
      ((mask[i] & groupbit2) && (mask[j] & groupbit))) {

    if (newton || i < nlocal) {
      if (mask[i] & groupbit) {
        ftotal[0] += fpair * dx;
        ftotal[1] += fpair * dy;
        ftotal[2] += fpair * dz;
      }
      fatom[i][0] += fpair * dx;
      fatom[i][1] += fpair * dy;
      fatom[i][2] += fpair * dz;
    }
    if (newton || j < nlocal) {
      if (mask[j] & groupbit) {
        ftotal[0] -= fpair * dx;
        ftotal[1] -= fpair * dy;
        ftotal[2] -= fpair * dz;
      }
      fatom[j][0] -= fpair * dx;
      fatom[j][1] -= fpair * dy;
      fatom[j][2] -= fpair * dz;
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeForceTally::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = fatom[i][0];
    buf[m++] = fatom[i][1];
    buf[m++] = fatom[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeForceTally::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    fatom[j][0] += buf[m++];
    fatom[j][1] += buf[m++];
    fatom[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double ComputeForceTally::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if ((did_setup != invoked_scalar) || (update->eflag_global != invoked_scalar))
    error->all(FLERR, "Energy was not tallied on needed timestep");

  // sum accumulated forces across procs

  MPI_Allreduce(ftotal, vector, size_peratom_cols, MPI_DOUBLE, MPI_SUM, world);

  scalar = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeForceTally::compute_peratom()
{
  invoked_peratom = update->ntimestep;
  if ((did_setup != invoked_peratom) || (update->eflag_global != invoked_peratom))
    error->all(FLERR, "Energy was not tallied on needed timestep");

  // collect contributions from ghost atoms

  if (force->newton_pair) {
    comm->reverse_comm(this);

    // clear out ghost atom data after it has been collected to local atoms
    const int nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; ++i)
      for (int j = 0; j < size_peratom_cols; ++j) fatom[i][j] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeForceTally::memory_usage()
{
  double bytes = (nmax < 0) ? 0 : nmax * (double) size_peratom_cols * sizeof(double);
  return bytes;
}

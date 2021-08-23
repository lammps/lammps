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

#include "compute_stress_tally.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeStressTally::ComputeStressTally(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR, "Illegal compute stress/tally command");

  igroup2 = group->find(arg[3]);
  if (igroup2 == -1) error->all(FLERR, "Could not find compute stress/tally second group ID");
  groupbit2 = group->bitmask[igroup2];

  scalar_flag = 1;
  vector_flag = 0;
  peratom_flag = 1;
  timeflag = 1;
  dynamic_group_allow = 0;

  comm_reverse = size_peratom_cols = 6;
  extscalar = 0;
  peflag = 1;    // we need Pair::ev_tally() to be run

  did_setup = invoked_peratom = invoked_scalar = -1;
  nmax = -1;
  stress = nullptr;
  vector = new double[size_peratom_cols];
  virial = new double[size_peratom_cols];
}

/* ---------------------------------------------------------------------- */

ComputeStressTally::~ComputeStressTally()
{
  if (force && force->pair) force->pair->del_tally_callback(this);
  memory->destroy(stress);
  delete[] virial;
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeStressTally::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Trying to use compute stress/tally without pair style");
  else
    force->pair->add_tally_callback(this);

  if (comm->me == 0) {
    if (force->pair->single_enable == 0 || force->pair->manybody_flag)
      error->warning(FLERR, "Compute stress/tally used with incompatible pair style");

    if (force->bond || force->angle || force->dihedral || force->improper || force->kspace)
      error->warning(FLERR, "Compute stress/tally only called from pair style");
  }
  did_setup = -1;
}

/* ---------------------------------------------------------------------- */

void ComputeStressTally::pair_setup_callback(int, int)
{
  // run setup only once per time step.
  // we may be called from multiple pair styles

  if (did_setup == update->ntimestep) return;

  const int ntotal = atom->nlocal + atom->nghost;

  // grow per-atom storage, if needed

  if (atom->nmax > nmax) {
    memory->destroy(stress);
    nmax = atom->nmax;
    memory->create(stress, nmax, size_peratom_cols, "stress/tally:stress");
    array_atom = stress;
  }

  // clear storage

  for (int i = 0; i < ntotal; ++i)
    for (int j = 0; j < size_peratom_cols; ++j) stress[i][j] = 0.0;

  for (int i = 0; i < size_peratom_cols; ++i) vector[i] = virial[i] = 0.0;

  did_setup = update->ntimestep;
}

/* ---------------------------------------------------------------------- */
void ComputeStressTally::pair_tally_callback(int i, int j, int nlocal, int newton, double, double,
                                             double fpair, double dx, double dy, double dz)
{
  const int *const mask = atom->mask;

  if (((mask[i] & groupbit) && (mask[j] & groupbit2)) ||
      ((mask[i] & groupbit2) && (mask[j] & groupbit))) {

    fpair *= 0.5;
    const double v0 = dx * dx * fpair;
    const double v1 = dy * dy * fpair;
    const double v2 = dz * dz * fpair;
    const double v3 = dx * dy * fpair;
    const double v4 = dx * dz * fpair;
    const double v5 = dy * dz * fpair;

    if (newton || i < nlocal) {
      virial[0] += v0;
      stress[i][0] += v0;
      virial[1] += v1;
      stress[i][1] += v1;
      virial[2] += v2;
      stress[i][2] += v2;
      virial[3] += v3;
      stress[i][3] += v3;
      virial[4] += v4;
      stress[i][4] += v4;
      virial[5] += v5;
      stress[i][5] += v5;
    }
    if (newton || j < nlocal) {
      virial[0] += v0;
      stress[j][0] += v0;
      virial[1] += v1;
      stress[j][1] += v1;
      virial[2] += v2;
      stress[j][2] += v2;
      virial[3] += v3;
      stress[j][3] += v3;
      virial[4] += v4;
      stress[j][4] += v4;
      virial[5] += v5;
      stress[j][5] += v5;
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeStressTally::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = stress[i][0];
    buf[m++] = stress[i][1];
    buf[m++] = stress[i][2];
    buf[m++] = stress[i][3];
    buf[m++] = stress[i][4];
    buf[m++] = stress[i][5];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeStressTally::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    stress[j][0] += buf[m++];
    stress[j][1] += buf[m++];
    stress[j][2] += buf[m++];
    stress[j][3] += buf[m++];
    stress[j][4] += buf[m++];
    stress[j][5] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double ComputeStressTally::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if ((did_setup != invoked_scalar) || (update->eflag_global != invoked_scalar))
    error->all(FLERR, "Energy was not tallied on needed timestep");

  // sum accumulated forces across procs

  MPI_Allreduce(virial, vector, size_peratom_cols, MPI_DOUBLE, MPI_SUM, world);

  if (domain->dimension == 3)
    scalar = (vector[0] + vector[1] + vector[2]) / 3.0;
  else
    scalar = (vector[0] + vector[1]) / 2.0;

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeStressTally::compute_peratom()
{
  invoked_peratom = update->ntimestep;
  if ((did_setup != invoked_peratom) || (update->eflag_global != invoked_peratom))
    error->all(FLERR, "Energy was not tallied on needed timestep");

  // collect contributions from ghost atoms

  if (force->newton_pair) {
    comm->reverse_comm_compute(this);

    const int nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; ++i)
      for (int j = 0; j < size_peratom_cols; ++j) stress[i][j] = 0.0;
  }

  // convert to stress*volume units = -pressure*volume

  const double nktv2p = -force->nktv2p;
  for (int i = 0; i < atom->nlocal; i++) {
    stress[i][0] *= nktv2p;
    stress[i][1] *= nktv2p;
    stress[i][2] *= nktv2p;
    stress[i][3] *= nktv2p;
    stress[i][4] *= nktv2p;
    stress[i][5] *= nktv2p;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeStressTally::memory_usage()
{
  double bytes = (nmax < 0) ? 0 : nmax * (double)size_peratom_cols * sizeof(double);
  return bytes;
}

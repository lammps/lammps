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

#include "compute_heat_flux_virial_tally.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeHeatFluxVirialTally::ComputeHeatFluxVirialTally(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR, "Illegal compute heat/flux/virial/tally command");

  igroup2 = group->find(arg[3]);
  if (igroup2 == -1)
    error->all(FLERR, "Could not find compute heat/flux/virial/tally second group ID");
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
}

/* ---------------------------------------------------------------------- */

ComputeHeatFluxVirialTally::~ComputeHeatFluxVirialTally()
{
  if (force && force->pair) force->pair->del_tally_callback(this);
  memory->destroy(fatom);
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVirialTally::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Trying to use compute heat/flux/virial/tally without pair style");
  else
    force->pair->add_tally_callback(this);

  if (comm->me == 0) {
    if (force->pair->single_enable == 0 || force->pair->manybody_flag)
      error->warning(FLERR, "Compute heat/flux/virial/tally used with incompatible pair style");

    if (force->bond || force->angle || force->dihedral || force->improper || force->kspace)
      error->warning(FLERR, "Compute heat/flux/virial/tally only called from pair style");
  }

  // error out if any atoms are in both groups
  for (int i = 0; i < atom->nlocal; i++) {
    if ((atom->mask[i] & groupbit) && (atom->mask[i] & groupbit2)) {
      if (atom->tag_enable) {
        error->all(FLERR,
                   "Atom {} belonging to both groups is not allowed"
                   " with compute heat/flux/virial/tally",
                   atom->tag[i]);
      } else {
        error->all(FLERR,
                   "Atom belonging to both groups is not allowed"
                   " with compute heat/flux/virial/tally");
      }
    }
  }

  did_setup = -1;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVirialTally::pair_setup_callback(int, int)
{
  // run setup only once per time step.
  // we may be called from multiple pair styles

  if (did_setup == update->ntimestep) return;

  const int ntotal = atom->nlocal + atom->nghost;

  // grow per-atom storage, if needed

  if (atom->nmax > nmax) {
    memory->destroy(fatom);
    nmax = atom->nmax;
    memory->create(fatom, nmax, size_peratom_cols, "heat/flux/virial/tally:fatom");
    array_atom = fatom;
  }

  // clear storage

  for (int i = 0; i < ntotal; ++i)
    for (int j = 0; j < size_peratom_cols; ++j) fatom[i][j] = 0.0;

  did_setup = update->ntimestep;
}

/* ---------------------------------------------------------------------- */
void ComputeHeatFluxVirialTally::pair_tally_callback(int i, int j, int nlocal, int newton, double,
                                                     double, double fpair, double dx, double dy,
                                                     double dz)
{
  const int *const mask = atom->mask;

  const bool ingroup1 = (mask[i] & groupbit);
  if ((ingroup1 && (mask[j] & groupbit2)) || ((mask[i] & groupbit2) && (mask[j] & groupbit))) {

    // signs set to calculate heat flux from group2 into group1
    const double fx = (ingroup1 ? 0.5 : -0.5) * dx * fpair;
    const double fy = (ingroup1 ? 0.5 : -0.5) * dy * fpair;
    const double fz = (ingroup1 ? 0.5 : -0.5) * dz * fpair;

    if (newton || i < nlocal) {
      fatom[i][0] += fx;
      fatom[i][1] += fy;
      fatom[i][2] += fz;
    }
    if (newton || j < nlocal) {
      fatom[j][0] += fx;
      fatom[j][1] += fy;
      fatom[j][2] += fz;
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeHeatFluxVirialTally::pack_reverse_comm(int n, int first, double *buf)
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

void ComputeHeatFluxVirialTally::unpack_reverse_comm(int n, int *list, double *buf)
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

double ComputeHeatFluxVirialTally::compute_scalar()
{
  // need to collect contributions from ghost atoms
  // because we need to use local velocities to compute heat flux
  if (invoked_peratom != update->ntimestep) compute_peratom();

  invoked_scalar = update->ntimestep;
  if ((did_setup != invoked_scalar) || (update->eflag_global != invoked_scalar))
    error->all(FLERR, "Energy was not tallied on needed timestep");

  // sum heat flux across procs
  double hflux = 0.0;
  for (int i = 0; i < atom->nlocal; i++)
    hflux +=
        fatom[i][0] * atom->v[i][0] + fatom[i][1] * atom->v[i][1] + fatom[i][2] * atom->v[i][2];

  MPI_Allreduce(&hflux, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFluxVirialTally::compute_peratom()
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

double ComputeHeatFluxVirialTally::memory_usage()
{
  double bytes = (nmax < 0) ? 0 : nmax * (double) size_peratom_cols * sizeof(double);
  return bytes;
}

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

#include "compute_pe_tally.h"

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

ComputePETally::ComputePETally(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "compute pe/tally", error);

  igroup2 = group->find(arg[3]);
  if (igroup2 == -1) error->all(FLERR, "Could not find compute pe/tally second group ID {}", arg[3]);
  groupbit2 = group->bitmask[igroup2];

  scalar_flag = 1;
  vector_flag = 0;
  peratom_flag = 1;
  timeflag = 1;

  comm_reverse = size_peratom_cols = 2;
  extscalar = 1;
  peflag = 1;    // we need Pair::ev_tally() to be run

  did_setup = invoked_peratom = invoked_scalar = -1;
  nmax = -1;
  eatom = nullptr;
  vector = new double[size_peratom_cols];
}

/* ---------------------------------------------------------------------- */

ComputePETally::~ComputePETally()
{
  if (force && force->pair) force->pair->del_tally_callback(this);
  memory->destroy(eatom);
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputePETally::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Trying to use compute pe/tally without a pair style");
  else
    force->pair->add_tally_callback(this);

  if (comm->me == 0) {
    if (force->pair->single_enable == 0 || force->pair->manybody_flag)
      error->warning(FLERR, "Compute pe/tally used with incompatible pair style");

    if (force->bond || force->angle || force->dihedral || force->improper || force->kspace)
      error->warning(FLERR, "Compute pe/tally only called from pair style");
  }
  did_setup = -1;
}

/* ---------------------------------------------------------------------- */

void ComputePETally::pair_setup_callback(int, int)
{
  // run setup only once per time step.
  // we may be called from multiple pair styles

  if (did_setup == update->ntimestep) return;

  const int ntotal = atom->nlocal + atom->nghost;

  // grow per-atom storage, if needed

  if (atom->nmax > nmax) {
    memory->destroy(eatom);
    nmax = atom->nmax;
    memory->create(eatom, nmax, size_peratom_cols, "pe/tally:eatom");
    array_atom = eatom;
  }

  // clear storage

  for (int i = 0; i < ntotal; ++i) eatom[i][0] = eatom[i][1] = 0.0;

  vector[0] = etotal[0] = vector[1] = etotal[1] = 0.0;

  did_setup = update->ntimestep;
}

/* ---------------------------------------------------------------------- */
void ComputePETally::pair_tally_callback(int i, int j, int nlocal, int newton, double evdwl,
                                         double ecoul, double, double, double, double)
{
  const int *const mask = atom->mask;

  if (((mask[i] & groupbit) && (mask[j] & groupbit2)) ||
      ((mask[i] & groupbit2) && (mask[j] & groupbit))) {

    evdwl *= 0.5;
    ecoul *= 0.5;
    if (newton || i < nlocal) {
      etotal[0] += evdwl;
      eatom[i][0] += evdwl;
      etotal[1] += ecoul;
      eatom[i][1] += ecoul;
    }
    if (newton || j < nlocal) {
      etotal[0] += evdwl;
      eatom[j][0] += evdwl;
      etotal[1] += ecoul;
      eatom[j][1] += ecoul;
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputePETally::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = eatom[i][0];
    buf[m++] = eatom[i][1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputePETally::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    eatom[j][0] += buf[m++];
    eatom[j][1] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double ComputePETally::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if ((did_setup != invoked_scalar) || (update->eflag_global != invoked_scalar))
    error->all(FLERR, "Energy was not tallied on needed timestep");

  if ((comm->me == 0) && !force->pair->did_tally_callback())
    error->warning(FLERR, "Energy was not tallied by pair style");

  // sum accumulated energies across procs

  MPI_Allreduce(etotal, vector, size_peratom_cols, MPI_DOUBLE, MPI_SUM, world);

  scalar = vector[0] + vector[1];
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputePETally::compute_peratom()
{
  invoked_peratom = update->ntimestep;
  if ((did_setup != invoked_peratom) || (update->eflag_global != invoked_peratom))
    error->all(FLERR, "Energy was not tallied on needed timestep");

  if ((comm->me == 0) && !force->pair->did_tally_callback())
    error->warning(FLERR, "Energy was not tallied by pair style");

  // collect contributions from ghost atoms

  if (force->newton_pair) {
    comm->reverse_comm(this);

    // clear out ghost atom data after it has been collected to local atoms
    const int nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; ++i) eatom[i][0] = eatom[i][1] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePETally::memory_usage()
{
  double bytes = (nmax < 0) ? 0 : nmax * (double) size_peratom_cols * sizeof(double);
  return bytes;
}

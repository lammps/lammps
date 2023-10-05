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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "compute_aggregate_atom.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeAggregateAtom::ComputeAggregateAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), aggregateID(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute aggregate/atom command");

  double cutoff = utils::numeric(FLERR, arg[3], false, lmp);
  cutsq = cutoff * cutoff;

  if (atom->avec->bonds_allow == 0)
    error->all(FLERR, "Compute aggregate/atom used when bonds are not allowed");

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_forward = 1;
  comm_reverse = 1;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeAggregateAtom::~ComputeAggregateAtom()
{
  memory->destroy(aggregateID);
}

/* ---------------------------------------------------------------------- */

void ComputeAggregateAtom::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR, "Cannot use compute aggregate/atom unless atoms have IDs");
  if (force->bond == nullptr)
    error->all(FLERR, "Compute aggregate/atom requires a bond style to be defined");

  if (force->pair == nullptr)
    error->all(FLERR, "Compute cluster/atom requires a pair style to be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR, "Compute cluster/atom cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list
  // full required so that pair of atoms on 2 procs both set their clusterID

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style(style).size() > 1)
    if (comm->me == 0) error->warning(FLERR, "More than one compute {}", style);
}

/* ---------------------------------------------------------------------- */

void ComputeAggregateAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeAggregateAtom::compute_peratom()
{
  int i, j, k;

  invoked_peratom = update->ntimestep;

  // grow aggregateID array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(aggregateID);
    nmax = atom->nmax;
    memory->create(aggregateID, nmax, "aggregate/atom:aggregateID");
    vector_atom = aggregateID;
  }

  // invoke full neighbor list (will copy or build if necessary)
  // on the first step of a run, set preflag to one in neighbor->build_one(...)

  if (update->firststep == update->ntimestep)
    neighbor->build_one(list, 1);
  else
    neighbor->build_one(list);

  // if group is dynamic, ensure ghost atom masks are current

  if (group->dynamic[igroup]) {
    commflag = 0;
    comm->forward_comm(this);
  }

  // each atom starts in its own aggregate,

  int nlocal = atom->nlocal;
  int inum = list->inum;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  double **x = atom->x;

  for (i = 0; i < nlocal + atom->nghost; i++)
    if (mask[i] & groupbit)
      aggregateID[i] = tag[i];
    else
      aggregateID[i] = 0;

  // loop until no more changes on any proc:
  // acquire aggregateIDs of ghost atoms
  // loop over my atoms, and check atoms bound to it
  // if both atoms are in aggregate, assign lowest aggregateID to both
  // then loop over my atoms, checking distance to neighbors
  // if both atoms are in cluster, assign lowest clusterID to both
  // iterate until no changes in my atoms
  // then check if any proc made changes

  commflag = 1;

  int change, done, anychange;

  while (true) {
    comm->forward_comm(this);

    // reverse communication when bonds are not stored on every processor

    if (force->newton_bond) comm->reverse_comm(this);

    change = 0;
    while (true) {
      done = 1;
      for (i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;

        for (j = 0; j < num_bond[i]; j++) {
          if (bond_type[i][j] == 0) continue;
          k = atom->map(bond_atom[i][j]);
          if (k < 0) continue;
          if (!(mask[k] & groupbit)) continue;
          if (aggregateID[i] == aggregateID[k]) continue;

          aggregateID[i] = aggregateID[k] = MIN(aggregateID[i], aggregateID[k]);
          done = 0;
        }
      }

      for (int ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        if (!(mask[i] & groupbit)) continue;

        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        int *jlist = firstneigh[i];
        const int jnum = numneigh[i];

        for (int jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          if (!(mask[j] & groupbit)) continue;
          if (aggregateID[i] == aggregateID[j]) continue;

          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            aggregateID[i] = aggregateID[j] = MIN(aggregateID[i], aggregateID[j]);
            done = 0;
          }
        }
      }
      if (!done) change = 1;
      if (done) break;
    }

    // stop if all procs are done

    MPI_Allreduce(&change, &anychange, 1, MPI_INT, MPI_MAX, world);
    if (!anychange) break;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeAggregateAtom::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                            int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  if (commflag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = aggregateID[j];
    }
  } else {
    int *mask = atom->mask;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(mask[j]).d;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeAggregateAtom::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  if (commflag)
    for (i = first; i < last; i++) {
      double x = buf[m++];

      // only overwrite ghost IDs with values lower than current ones

      aggregateID[i] = MIN(x, aggregateID[i]);
    }
  else {
    int *mask = atom->mask;
    for (i = first; i < last; i++) mask[i] = (int) ubuf(buf[m++]).i;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeAggregateAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) { buf[m++] = aggregateID[i]; }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeAggregateAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    double x = buf[m++];

    // only overwrite local IDs with values lower than current ones

    aggregateID[j] = MIN(x, aggregateID[j]);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeAggregateAtom::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}

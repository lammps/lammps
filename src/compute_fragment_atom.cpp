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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include <cstring>
#include "compute_fragment_atom.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeFragmentAtom::ComputeFragmentAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  fragmentID(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute fragment/atom command");

  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Compute fragment/atom used when bonds are not allowed");

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_forward = 1;
  comm_reverse = 1;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeFragmentAtom::~ComputeFragmentAtom()
{
  memory->destroy(fragmentID);
}

/* ---------------------------------------------------------------------- */

void ComputeFragmentAtom::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute fragment/atom unless atoms have IDs");
  if (force->bond == NULL)
    error->all(FLERR,"Compute fragment/atom requires a bond style to be defined");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"fragment/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute fragment/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeFragmentAtom::compute_peratom()
{
  int i,j,k;

  invoked_peratom = update->ntimestep;

  // grow fragmentID array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(fragmentID);
    nmax = atom->nmax;
    memory->create(fragmentID,nmax,"fragment/atom:fragmentID");
    vector_atom = fragmentID;
  }

  // if group is dynamic, insure ghost atom masks are current

  if (group->dynamic[igroup]) {
    commflag = 0;
    comm->forward_comm_compute(this);
  }

  // each atom starts in its own fragment,

  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;

  for (i = 0; i < nlocal + atom->nghost; i++)
    if (mask[i] & groupbit) fragmentID[i] = tag[i];
    else fragmentID[i] = 0;

  // loop until no more changes on any proc:
  // acquire fragmentIDs of ghost atoms
  // loop over my atoms, and check atoms bound to it
  // if both atoms are in fragment, assign lowest fragmentID to both
  // iterate until no changes in my atoms
  // then check if any proc made changes

  commflag = 1;

  int change,done,anychange;

  while (1) {
    comm->forward_comm_compute(this);

    // reverse communication when bonds are not stored on every processor

    if (force->newton_bond)
      comm->reverse_comm_compute(this);

    change = 0;
    while (1) {
      done = 1;
      for (i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;

        for (j = 0; j < num_bond[i]; j++) {
          if (bond_type[i][j] == 0) continue;
          k = atom->map(bond_atom[i][j]);
          if (k < 0) continue;
          if (!(mask[k] & groupbit)) continue;
          if (fragmentID[i] == fragmentID[k]) continue;

          fragmentID[i] = fragmentID[k] = MIN(fragmentID[i],fragmentID[k]);
          done = 0;
        }
      }
      if (!done) change = 1;
      if (done) break;
    }

    // stop if all procs are done

    MPI_Allreduce(&change,&anychange,1,MPI_INT,MPI_MAX,world);
    if (!anychange) break;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeFragmentAtom::pack_forward_comm(int n, int *list, double *buf,
                                          int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  if (commflag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = fragmentID[j];
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

void ComputeFragmentAtom::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag)
    for (i = first; i < last; i++) {
      double x = buf[m++];

      // only overwrite ghost IDs with values lower than current ones

      fragmentID[i] = MIN(x,fragmentID[i]);
    }
  else {
    int *mask = atom->mask;
    for (i = first; i < last; i++) mask[i] = (int) ubuf(buf[m++]).i;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeFragmentAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = fragmentID[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeFragmentAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    double x = buf[m++];

    // only overwrite local IDs with values lower than current ones

    fragmentID[j] = MIN(x,fragmentID[j]);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeFragmentAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}

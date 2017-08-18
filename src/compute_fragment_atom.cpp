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

#include <string.h>
#include "compute_fragment_atom.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
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

  int nbondlist = neighbor->nbondlist;
  int **bondlist = neighbor->bondlist;

  // if group is dynamic, insure ghost atom masks are current

  if (group->dynamic[igroup]) {
    commflag = 0;
    comm->forward_comm_compute(this);
  }

  // every bond starts in its own fragment,
  // with fragmentID = MIN(b1atomID,b2atomID)
  // only bonds wholly contained in the group are considered

  tagint *tag = atom->tag;
  int *mask = atom->mask;

  for (i = 0; i < nbondlist; i++) {
    const int b1 = bondlist[i][0];
    const int b2 = bondlist[i][1];

    if ((mask[b1] & groupbit) && (mask[b2] & groupbit))
      fragmentID[b1] = fragmentID[b2] = MIN(tag[b1],tag[b2]);
    else fragmentID[b1] = fragmentID[b2] = 0;
  }

  // loop until no more changes on any proc:
  // acquire fragmentIDs of ghost atoms
  // loop over my atoms, and check atoms bound to it
  // if both atoms are in fragment, assign lowest fragmentID to both
  // iterate until no changes in my atoms
  // then check if any proc made changes

  commflag = 1;
  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;

  int change,done,anychange;

  while (1) {
    comm->forward_comm_compute(this);

    change = 0;
    while (1) {
      done = 1;
      for (i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        
        for (j = 0; j < num_bond[i]; j++) {
          k = bond_atom[i][j];
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
                                          int pbc_flag, int *pbc)
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
    for (i = first; i < last; i++) fragmentID[i] = buf[m++];
  else {
    int *mask = atom->mask;
    for (i = first; i < last; i++) mask[i] = (int) ubuf(buf[m++]).i;
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

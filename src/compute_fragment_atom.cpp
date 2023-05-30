// clang-format off
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

#include "compute_fragment_atom.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeFragmentAtom::ComputeFragmentAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  fragmentID(nullptr)
{
  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Compute fragment/atom used when bonds are not allowed");

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_forward = 1;

  // process optional args

  singleflag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"single") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute fragment/atom command");
      singleflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute fragment/atom command");
  }

  nmax = 0;
  stack = nullptr;
  clist = nullptr;
  markflag = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeFragmentAtom::~ComputeFragmentAtom()
{
  memory->destroy(stack);
  memory->destroy(clist);
  memory->destroy(markflag);
  memory->destroy(fragmentID);
}

/* ---------------------------------------------------------------------- */

void ComputeFragmentAtom::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute fragment/atom unless atoms have IDs");
  if (atom->molecular != Atom::MOLECULAR)
    error->all(FLERR,"Compute fragment/atom requires a molecular system");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"fragment/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute fragment/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeFragmentAtom::compute_peratom()
{
  int i,j,k,m,n;
  int nstack,ncluster,done,alldone;
  double newID,cID;
  tagint *list;

  invoked_peratom = update->ntimestep;

  // grow work and fragmentID vectors if necessary

  if (atom->nmax > nmax) {
    memory->destroy(stack);
    memory->destroy(clist);
    memory->destroy(markflag);
    memory->destroy(fragmentID);
    nmax = atom->nmax;
    memory->create(stack,nmax,"fragment/atom:stack");
    memory->create(clist,nmax,"fragment/atom:clist");
    memory->create(markflag,nmax,"fragment/atom:markflag");
    memory->create(fragmentID,nmax,"fragment/atom:fragmentID");
    vector_atom = fragmentID;
  }

  // if group is dynamic, ensure ghost atom masks are current

  if (group->dynamic[igroup]) {
    commflag = 0;
    comm->forward_comm(this);
  }

  // owned + ghost atoms start with fragmentID = atomID
  // atoms not in group have fragmentID = 0

  tagint *tag = atom->tag;
  int *mask = atom->mask;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    if (mask[i] & groupbit) fragmentID[i] = tag[i];
    else fragmentID[i] = 0;
  }

  // loop until no ghost atom fragment ID is changed
  // acquire fragmentIDs of ghost atoms
  // loop over clusters of atoms, which include ghost atoms
  // set fragmentIDs for each cluster to min framentID in the clusters
  // if singleflag = 0 atoms without bonds are assigned fragmentID = 0
  // iterate until no changes to ghost atom fragmentIDs

  commflag = 1;

  while (true) {
    comm->forward_comm(this);
    done = 1;

    // set markflag = 0 for all owned atoms, for new iteration

    for (i = 0; i < nlocal; i++) markflag[i] = 0;

    // loop over all owned atoms
    // each unmarked atom starts a cluster search

    for (i = 0; i < nlocal; i++) {

      // skip atom I if not in group or already marked
      // if singleflag = 0 and atom has no bond partners, fragID = 0 and done

      if (!(mask[i] & groupbit)) continue;
      if (markflag[i]) continue;
      if (!singleflag && (nspecial[i][0] == 0)) {
        fragmentID[i] = 0.0;
        continue;
      }

      // find one cluster of bond-connected atoms
      // ncluster = # of owned and ghost atoms in cluster
      // clist = vector of local indices of the ncluster atoms
      // stack is used to walk the bond topology

      ncluster = nstack = 0;
      stack[nstack++] = i;

      while (nstack) {
        j = stack[--nstack];
        clist[ncluster++] = j;
        markflag[j] = 1;

        n = nspecial[j][0];
        list = special[j];
        for (m = 0; m < n; m++) {
          k = atom->map(list[m]);

          // skip bond neighbor K if not in group or already marked

          if (k < 0) continue;
          if (!(mask[k] & groupbit)) continue;
          if (k < nlocal && markflag[k]) continue;

          // owned bond neighbors are added to stack for further walking
          // ghost bond neighbors are added directly w/out use of stack

          if (k < nlocal) stack[nstack++] = k;
          else clist[ncluster++] = k;
        }
      }

      // newID = minimum fragment ID in cluster list, including ghost atoms

      newID = BIG;
      for (m = 0; m < ncluster; m++) {
        cID = fragmentID[clist[m]];
        newID = MIN(newID,cID);
      }

      // set fragmentID = newID for all atoms in cluster, including ghost atoms
      // not done with iterations if change the fragmentID of a ghost atom

      for (m = 0; m < ncluster; m++) {
        j = clist[m];
        if (j >= nlocal && fragmentID[j] != newID) done = 0;
        fragmentID[j] = newID;
      }
    }

    // stop if all procs are done

    MPI_Allreduce(&done,&alldone,1,MPI_INT,MPI_MIN,world);
    if (alldone) break;
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

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double ComputeFragmentAtom::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  bytes += (double)3*nmax * sizeof(int);
  return bytes;
}

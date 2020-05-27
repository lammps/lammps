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

#include "npair_skip.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairSkip::NPairSkip(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   works for half and full lists
   works for owned (non-ghost) list, also for ghost list
   iskip and ijskip flag which atom types and type pairs to skip
   if ghost, also store neighbors of ghost atoms & set inum,gnum correctly
------------------------------------------------------------------------- */

void NPairSkip::build(NeighList *list)
{
  int i,j,ii,jj,n,itype,jnum,joriginal;
  int *neighptr,*jlist;

  int *type = atom->type;
  int nlocal = atom->nlocal;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int *ilist_skip = list->listskip->ilist;
  int *numneigh_skip = list->listskip->numneigh;
  int **firstneigh_skip = list->listskip->firstneigh;
  int num_skip = list->listskip->inum;
  if (list->ghost) num_skip += list->listskip->gnum;

  int *iskip = list->iskip;
  int **ijskip = list->ijskip;

  int inum = 0;
  ipage->reset();

  // loop over atoms in other list
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (ii = 0; ii < num_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[itype]) continue;

    n = 0;
    neighptr = ipage->vget();

    // loop over parent non-skip list

    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (ijskip[itype][type[j]]) continue;
      neighptr[n++] = joriginal;
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
  if (list->ghost) {
    int num = 0;
    for (i = 0; i < inum; i++)
      if (ilist[i] < nlocal) num++;
      else break;
    list->inum = num;
    list->gnum = inum - num;
  }
}

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
   Contributing authors:
   James Larentzos and Timothy I. Mattox (Engility Corporation)
------------------------------------------------------------------------- */

#include "npair_halffull_newton_ssa.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

// allocate space for static class variable
// prototype for non-class function

static int *ssaAIRptr;
static int cmp_ssaAIR(const void *, const void *);

/* ---------------------------------------------------------------------- */

NPairHalffullNewtonSSA::NPairHalffullNewtonSSA(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   build half list from full list for use by Shardlow Spliting Algorithm
   pair stored once if i,j are both owned and i < j
   if j is ghost, only store if j coords are "above and to the right" of i
   works if full list is a skip list
------------------------------------------------------------------------- */

void NPairHalffullNewtonSSA::build(NeighList *list)
{
  int i,j,ii,jj,n,jnum,joriginal;
  int *neighptr,*jlist;

  int nlocal = atom->nlocal;
  int *ssaAIR = atom->ssaAIR;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int *ilist_full = list->listfull->ilist;
  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;
  int inum_full = list->listfull->inum;

  int inum = 0;
  ipage->reset();

  // loop over parent full list

  for (ii = 0; ii < inum_full; ii++) {
    int AIRct[8] = { 0 };
    n = 0;
    neighptr = ipage->vget();

    i = ilist_full[ii];

    // loop over full neighbor list

    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (j < nlocal) {
        if (i > j) continue;
        ++(AIRct[0]);
      } else {
        if (ssaAIR[j] < 2) continue; // skip ghost atoms not in AIR
        ++(AIRct[ssaAIR[j] - 1]);
      }
      neighptr[n++] = joriginal;
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    // sort the locals+ghosts in the neighbor list by their ssaAIR number

    ssaAIRptr = atom->ssaAIR;
    qsort(&(neighptr[0]), n, sizeof(int), cmp_ssaAIR);

    // do a prefix sum on the counts to turn them into indexes

    list->ndxAIR_ssa[i][0] = AIRct[0];
    for (int ndx = 1; ndx < 8; ++ndx) {
      list->ndxAIR_ssa[i][ndx] = AIRct[ndx] + list->ndxAIR_ssa[i][ndx - 1];
    }
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   comparison function invoked by qsort()
   accesses static class member ssaAIRptr, set before call to qsort()
------------------------------------------------------------------------- */

static int cmp_ssaAIR(const void *iptr, const void *jptr)
{
  int i = NEIGHMASK & *((int *) iptr);
  int j = NEIGHMASK & *((int *) jptr);
  if (ssaAIRptr[i] < ssaAIRptr[j]) return -1;
  if (ssaAIRptr[i] > ssaAIRptr[j]) return 1;
  return 0;
}


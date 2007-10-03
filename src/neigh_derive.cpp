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

#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
   works if full list is a skip list
------------------------------------------------------------------------- */

void Neighbor::half_full_no_newton(NeighList *list)
{
  int i,j,ii,jj,n,jnum;
  int *neighptr,*jlist;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;
  int *ilist_full = list->listfull->ilist;
  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;
  int inum_full = list->listfull->inum;

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  // loop over atoms in full list

  for (ii = 0; ii < inum_full; ii++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }

    neighptr = &pages[npage][npnt];
    n = 0;

    // loop over full neighbor list

    i = ilist_full[ii];
    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];

    for (j = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (j > i) neighptr[n++] = j;
    }

    ilist[inum] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    inum++;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
   works if full list is a skip list
------------------------------------------------------------------------- */

void Neighbor::half_full_newton(NeighList *list)
{
  int i,j,ii,jj,n,jnum;
  int *neighptr,*jlist;
  double xtmp,ytmp,ztmp;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;
  int *ilist_full = list->listfull->ilist;
  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;
  int inum_full = list->listfull->inum;

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  // loop over atoms in full list

  for (ii = 0; ii < inum_full; ii++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }

    neighptr = &pages[npage][npnt];
    n = 0;

    i = ilist_full[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over full neighbor list

    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (j < nlocal) {
	if (i > j) continue;
      } else {
	if (x[j][2] < ztmp) continue;
	if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }
      neighptr[n++] = j;
    }

    ilist[inum] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    inum++;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   iskip and ijskip flag which atom types and type pairs to skip
------------------------------------------------------------------------- */

void Neighbor::skip_from(NeighList *list)
{
  int i,j,ii,jj,n,itype,jnum,joriginal;
  int *neighptr,*jlist;

  int *type = atom->type;
  int nall = atom->nlocal + atom->nghost;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;
  int *ilist_skip = list->listskip->ilist;
  int *numneigh_skip = list->listskip->numneigh;
  int **firstneigh_skip = list->listskip->firstneigh;
  int inum_skip = list->listskip->inum;

  int *iskip = list->iskip;
  int **ijskip = list->ijskip;

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  // loop over atoms in other list
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (ii = 0; ii < inum_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[type[i]]) continue;

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }

    neighptr = &pages[npage][npnt];
    n = 0;

    // loop over full neighbor list

    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      j = joriginal = jlist[jj];
      if (j >= nall) j %= nall;
      if (ijskip[itype][type[j]]) continue;
      neighptr[n++] = joriginal;
    }

    ilist[inum] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    inum++;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   create list which is simply a copy of another list
------------------------------------------------------------------------- */

void Neighbor::copy_from(NeighList *list)
{
  NeighList *listcopy = list->listcopy;

  list->inum = listcopy->inum;
  list->ilist = listcopy->ilist;
  list->numneigh = listcopy->numneigh;
  list->firstneigh = listcopy->firstneigh;
  list->firstdouble = listcopy->firstdouble;
  list->pages = listcopy->pages;
  list->dpages = listcopy->dpages;
}

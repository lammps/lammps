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

void Neighbor::half_from_full_no_newton(NeighList *list)
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

    // loop over parent full list

    i = ilist_full[ii];
    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (j > i) neighptr[n++] = j;
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom || npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i,j are both owned and i < j
   if j is ghost, only store if j coords are "above and to the right" of i
   works if full list is a skip list
------------------------------------------------------------------------- */

void Neighbor::half_from_full_newton(NeighList *list)
{
  int i,j,ii,jj,n,jnum,joriginal;
  int *neighptr,*jlist;
  double xtmp,ytmp,ztmp;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

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

  // loop over parent full list

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
      j = joriginal = jlist[jj];
      if (j < nlocal) {
	if (i > j) continue;
      } else {
	if (j >= nall) j %= nall;
	if (x[j][2] < ztmp) continue;
	if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }
      neighptr[n++] = joriginal;
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom || npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   iskip and ijskip flag which atom types and type pairs to skip
   this is for half and full lists
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
    if (iskip[itype]) continue;

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }

    neighptr = &pages[npage][npnt];
    n = 0;

    // loop over parent non-skip list

    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      j = joriginal = jlist[jj];
      if (j >= nall) j %= nall;
      if (ijskip[itype][type[j]]) continue;
      neighptr[n++] = joriginal;
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom || npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   iskip and ijskip flag which atom types and type pairs to skip
   this is for granular lists with history, copy the history values from parent
------------------------------------------------------------------------- */

void Neighbor::skip_from_granular(NeighList *list)
{
  int i,j,ii,jj,n,nn,itype,jnum,joriginal;
  int *neighptr,*jlist,*touchptr,*touchptr_skip;
  double *shearptr,*shearptr_skip;

  int *type = atom->type;
  int nall = atom->nlocal + atom->nghost;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;
  int *ilist_skip = list->listskip->ilist;
  int *numneigh_skip = list->listskip->numneigh;
  int **firstneigh_skip = list->listskip->firstneigh;
  int **firsttouch_skip = list->listskip->listgranhistory->firstneigh;
  double **firstshear_skip = list->listskip->listgranhistory->firstdouble;
  int inum_skip = list->listskip->inum;

  int *iskip = list->iskip;
  int **ijskip = list->ijskip;

  NeighList *listgranhistory = list->listgranhistory;
  int **firsttouch = listgranhistory->firstneigh;
  double **firstshear = listgranhistory->firstdouble;
  int **pages_touch = listgranhistory->pages;
  double **pages_shear = listgranhistory->dpages;

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  // loop over atoms in other list
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (ii = 0; ii < inum_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[itype]) continue;

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) {
	pages = list->add_pages();
	pages_touch = listgranhistory->add_pages();
	pages_shear = listgranhistory->dpages;
      }
    }

    n = 0;
    neighptr = &pages[npage][npnt];
    nn = 0;
    touchptr = &pages_touch[npage][npnt];
    shearptr = &pages_shear[npage][3*npnt];

    // loop over parent non-skip granular list and its history info

    touchptr_skip = firsttouch_skip[i];
    shearptr_skip = firstshear_skip[i];
    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      j = joriginal = jlist[jj];
      if (j >= nall) j %= nall;
      if (ijskip[itype][type[j]]) continue;
      neighptr[n] = joriginal;
      touchptr[n++] = touchptr_skip[jj];
      shearptr[nn++] = shearptr_skip[3*jj];
      shearptr[nn++] = shearptr_skip[3*jj+1];
      shearptr[nn++] = shearptr_skip[3*jj+2];
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    firsttouch[i] = touchptr;
    firstshear[i] = shearptr;
    npnt += n;
    if (n > oneatom || npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   iskip and ijskip flag which atom types and type pairs to skip
   this is for respa lists, copy the inner/middle values from parent
------------------------------------------------------------------------- */

void Neighbor::skip_from_respa(NeighList *list)
{
  int i,j,ii,jj,n,itype,jnum,joriginal,n_inner,n_middle;
  int *neighptr,*jlist,*neighptr_inner,*neighptr_middle;

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

  NeighList *listinner = list->listinner;
  int *numneigh_inner = listinner->numneigh;
  int **firstneigh_inner = listinner->firstneigh;
  int **pages_inner = listinner->pages;
  int *numneigh_inner_skip = list->listskip->listinner->numneigh;
  int **firstneigh_inner_skip = list->listskip->listinner->firstneigh;

  NeighList *listmiddle;
  int *numneigh_middle,**firstneigh_middle,**pages_middle;
  int *numneigh_middle_skip,**firstneigh_middle_skip;
  int respamiddle = list->respamiddle;
  if (respamiddle) {
    listmiddle = list->listmiddle;
    numneigh_middle = listmiddle->numneigh;
    firstneigh_middle = listmiddle->firstneigh;
    pages_middle = listmiddle->pages;
    numneigh_middle_skip = list->listskip->listmiddle->numneigh;
    firstneigh_middle_skip = list->listskip->listmiddle->firstneigh;
  }

  int inum = 0;
  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  // loop over atoms in other list
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (ii = 0; ii < inum_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[itype]) continue;

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }
    neighptr = &pages[npage][npnt];
    n = 0;

    if (pgsize - npnt_inner < oneatom) {
      npnt_inner = 0;
      npage_inner++;
      if (npage_inner == listinner->maxpage)
	pages_inner = listinner->add_pages();
    }
    neighptr_inner = &pages_inner[npage_inner][npnt_inner];
    n_inner = 0;

    if (respamiddle) {
      if (pgsize - npnt_middle < oneatom) {
	npnt_middle = 0;
	npage_middle++;
	if (npage_middle == listmiddle->maxpage)
	  pages_middle = listmiddle->add_pages();
      }
      neighptr_middle = &pages_middle[npage_middle][npnt_middle];
      n_middle = 0;
    }

    // loop over parent outer rRESPA list

    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      j = joriginal = jlist[jj];
      if (j >= nall) j %= nall;
      if (ijskip[itype][type[j]]) continue;
      neighptr[n++] = joriginal;
    }

    // loop over parent inner rRESPA list

    jlist = firstneigh_inner_skip[i];
    jnum = numneigh_inner_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      j = joriginal = jlist[jj];
      if (j >= nall) j %= nall;
      if (ijskip[itype][type[j]]) continue;
      neighptr_inner[n_inner++] = joriginal;
    }

    // loop over parent middle rRESPA list

    if (respamiddle) {
      jlist = firstneigh_middle_skip[i];
      jnum = numneigh_middle_skip[i];

      for (jj = 0; jj < jnum; jj++) {
	j = joriginal = jlist[jj];
	if (j >= nall) j %= nall;
	if (ijskip[itype][type[j]]) continue;
	neighptr_middle[n_middle++] = joriginal;
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom || npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (npnt_inner >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    if (respamiddle) {
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (npnt_middle >= pgsize)
	error->one("Neighbor list overflow, boost neigh_modify one or page");
    }
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   create list which is simply a copy of parent list
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

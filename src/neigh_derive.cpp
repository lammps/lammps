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
#include "my_page.h"
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
  int i,j,ii,jj,n,jnum,joriginal;
  int *neighptr,*jlist;

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

  // loop over atoms in full list

  for (ii = 0; ii < inum_full; ii++) {
    n = 0;
    neighptr = ipage->vget();

    // loop over parent full list

    i = ilist_full[ii];
    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (j > i) neighptr[n++] = joriginal;
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
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
    n = 0;
    neighptr = ipage->vget();

    i = ilist_full[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over full neighbor list

    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (j < nlocal) {
        if (i > j) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp) {
          if (x[j][1] < ytmp) continue;
          if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
        }
      }
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
}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   iskip and ijskip flag which atom types and type pairs to skip
   this is for half and full lists
   if ghostflag, also store neighbors of ghost atoms & set inum,gnum correctly
------------------------------------------------------------------------- */

void Neighbor::skip_from(NeighList *list)
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
  if (list->ghostflag) num_skip += list->listskip->gnum;

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
  if (list->ghostflag) {
    int num = 0;
    for (i = 0; i < inum; i++)
      if (ilist[i] < nlocal) num++;
      else break;
    list->inum = num;
    list->gnum = inum - num;
  }
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

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

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
  MyPage<int> *ipage_touch = listgranhistory->ipage;
  MyPage<double> *dpage_shear = listgranhistory->dpage;

  int inum = 0;
  ipage->reset();
  ipage_touch->reset();
  dpage_shear->reset();

  // loop over atoms in other list
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (ii = 0; ii < inum_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[itype]) continue;

    n = nn = 0;
    neighptr = ipage->vget();
    touchptr = ipage_touch->vget();
    shearptr = dpage_shear->vget();

    // loop over parent non-skip granular list and its history info

    touchptr_skip = firsttouch_skip[i];
    shearptr_skip = firstshear_skip[i];
    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
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
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    firsttouch[i] = touchptr;
    firstshear[i] = shearptr;
    ipage_touch->vgot(n);
    dpage_shear->vgot(nn);
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

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int *ilist_skip = list->listskip->ilist;
  int *numneigh_skip = list->listskip->numneigh;
  int **firstneigh_skip = list->listskip->firstneigh;
  int inum_skip = list->listskip->inum;

  int *iskip = list->iskip;
  int **ijskip = list->ijskip;

  NeighList *listinner = list->listinner;
  int *ilist_inner = listinner->ilist;
  int *numneigh_inner = listinner->numneigh;
  int **firstneigh_inner = listinner->firstneigh;
  MyPage<int> *ipage_inner = listinner->ipage;

  int *numneigh_inner_skip = list->listskip->listinner->numneigh;
  int **firstneigh_inner_skip = list->listskip->listinner->firstneigh;

  NeighList *listmiddle;
  int *ilist_middle,*numneigh_middle,**firstneigh_middle;
  MyPage<int> *ipage_middle;
  int *numneigh_middle_skip,**firstneigh_middle_skip;
  int respamiddle = list->respamiddle;
  if (respamiddle) {
    listmiddle = list->listmiddle;
    ilist_middle = listmiddle->ilist;
    numneigh_middle = listmiddle->numneigh;
    firstneigh_middle = listmiddle->firstneigh;
    ipage_middle = listmiddle->ipage;
    numneigh_middle_skip = list->listskip->listmiddle->numneigh;
    firstneigh_middle_skip = list->listskip->listmiddle->firstneigh;
  }

  int inum = 0;
  ipage->reset();
  ipage_inner->reset();
  if (respamiddle) ipage_middle->reset();

  // loop over atoms in other list
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (ii = 0; ii < inum_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[itype]) continue;

    n = n_inner = 0;
    neighptr = ipage->vget();
    neighptr_inner = ipage_inner->vget();
    if (respamiddle) {
      n_middle = 0;
      neighptr_middle = ipage_middle->vget();
    }

    // loop over parent outer rRESPA list

    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (ijskip[itype][type[j]]) continue;
      neighptr[n++] = joriginal;
    }

    // loop over parent inner rRESPA list

    jlist = firstneigh_inner_skip[i];
    jnum = numneigh_inner_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (ijskip[itype][type[j]]) continue;
      neighptr_inner[n_inner++] = joriginal;
    }

    // loop over parent middle rRESPA list

    if (respamiddle) {
      jlist = firstneigh_middle_skip[i];
      jnum = numneigh_middle_skip[i];

      for (jj = 0; jj < jnum; jj++) {
        joriginal = jlist[jj];
        j = joriginal & NEIGHMASK;
        if (ijskip[itype][type[j]]) continue;
        neighptr_middle[n_middle++] = joriginal;
      }
    }

    ilist[inum] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    ilist_inner[inum] = i;
    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    ipage_inner->vgot(n);
    if (ipage_inner->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    if (respamiddle) {
      ilist_middle[inum] = i;
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      ipage_middle->vgot(n);
      if (ipage_middle->status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    }

    inum++;
  }

  list->inum = inum;
  listinner->inum = inum;
  if (respamiddle) listmiddle->inum = inum;
}

/* ----------------------------------------------------------------------
   create list which is simply a copy of parent list
------------------------------------------------------------------------- */

void Neighbor::copy_from(NeighList *list)
{
  NeighList *listcopy = list->listcopy;

  list->inum = listcopy->inum;
  list->gnum = listcopy->gnum;
  list->ilist = listcopy->ilist;
  list->numneigh = listcopy->numneigh;
  list->firstneigh = listcopy->firstneigh;
  list->firstdouble = listcopy->firstdouble;
  list->ipage = listcopy->ipage;
  list->dpage = listcopy->dpage;
}

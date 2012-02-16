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
#include "neighbor_omp.h"
#include "neigh_list.h"
#include "atom.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
   works if full list is a skip list
------------------------------------------------------------------------- */

void Neighbor::half_from_full_no_newton_omp(NeighList *list)
{
  const int inum_full = list->listfull->inum;

  NEIGH_OMP_INIT;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NEIGH_OMP_SETUP(inum_full);

  int i,j,ii,jj,n,jnum,joriginal;
  int *neighptr,*jlist;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *ilist_full = list->listfull->ilist;
  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;

  // each thread works on its own page
  int npage = tid;
  int npnt = 0;

  // loop over atoms in full list

  for (ii = ifrom; ii < ito; ii++) {

#if defined(_OPENMP)
#pragma omp critical
#endif
    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage += nthreads;
      // only one thread at a time may check whether we 
      // need new neighbor list pages and then add to them.
      if (npage >= list->maxpage) list->add_pages(nthreads);
    }

    neighptr = &(list->pages[npage][npnt]);
    n = 0;

    // loop over parent full list

    i = ilist_full[ii];
    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (j > i) neighptr[n++] = joriginal;
    }

    ilist[ii] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
  NEIGH_OMP_CLOSE;
  list->inum = inum_full;
}

/* ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i,j are both owned and i < j
   if j is ghost, only store if j coords are "above and to the right" of i
   works if full list is a skip list
------------------------------------------------------------------------- */

void Neighbor::half_from_full_newton_omp(NeighList *list)
{
  const int inum_full = list->listfull->inum;

  NEIGH_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NEIGH_OMP_SETUP(inum_full);

  int i,j,ii,jj,n,jnum,joriginal;
  int *neighptr,*jlist;
  double xtmp,ytmp,ztmp;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *ilist_full = list->listfull->ilist;
  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;

  // each thread works on its own page
  int npage = tid;
  int npnt = 0;

  // loop over parent full list

  for (ii = ifrom; ii < ito; ii++) {

#if defined(_OPENMP)
#pragma omp critical
#endif
    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage += nthreads;
      // only one thread at a time may check  whether we
      // need new neighbor list pages and then add to them.
      if (npage >= list->maxpage) list->add_pages(nthreads);
    }

    neighptr = &(list->pages[npage][npnt]);
    n = 0;

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

    ilist[ii] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
  NEIGH_OMP_CLOSE;
  list->inum = inum_full;
}


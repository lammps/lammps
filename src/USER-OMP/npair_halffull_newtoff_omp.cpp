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

#include "npair_halffull_newtoff_omp.h"
#include "npair_omp.h"
#include "neigh_list.h"
#include "atom_vec.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalffullNewtoffOmp::NPairHalffullNewtoffOmp(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
   works if full list is a skip list
------------------------------------------------------------------------- */

void NPairHalffullNewtoffOmp::build(NeighList *list)
{
  const int inum_full = list->listfull->inum;

  NPAIR_OMP_INIT;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NPAIR_OMP_SETUP(inum_full);

  int i,j,ii,jj,n,jnum,joriginal;
  int *neighptr,*jlist;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *ilist_full = list->listfull->ilist;
  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  ipage.reset();

  // loop over atoms in full list

  for (ii = ifrom; ii < ito; ii++) {

    n = 0;
    neighptr = ipage.vget();

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
    ipage.vgot(n);
    if (ipage.status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
  NPAIR_OMP_CLOSE;
  list->inum = inum_full;
}

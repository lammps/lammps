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
#include "domain.h"
#include "group.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::half_nsq_no_newton_omp(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  const int bitmask = (includegroup) ? group->bitmask[includegroup] : 0;
  const int nall = atom->nlocal + atom->nghost;

  NEIGH_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NEIGH_OMP_SETUP(nlocal);

  int i,j,n,itype,jtype,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  int **special = atom->special;
  int **nspecial = atom->nspecial;
  int *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int molecular = atom->molecular;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  ipage.reset();

  // loop over owned atoms, storing neighbors

  for (i = ifrom; i < ito; i++) {

    n = 0;
    neighptr = ipage.vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over remaining atoms, owned and ghost
    // only store pair if i < j

    for (j = i+1; j < nall; j++) {
      if (includegroup && !(mask[j] & bitmask)) continue;
      jtype = type[j];
      if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq[itype][jtype]) {
        if (molecular) {
          which = find_special(special[i],nspecial[i],tag[j]);
          if (which == 0) neighptr[n++] = j;
          else if (domain->minimum_image_check(delx,dely,delz))
            neighptr[n++] = j;
          else if (which > 0) neighptr[n++] = j ^ (which << SBBITS);
        } else neighptr[n++] = j;
      }
    }

    ilist[i] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage.vgot(n);
    if (ipage.status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
  NEIGH_OMP_CLOSE;
  list->inum = nlocal;
}

/* ----------------------------------------------------------------------
   N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
   include neighbors of ghost atoms, but no "special neighbors" for ghosts
   pair stored once if i,j are both owned and i < j
   pair stored by me if i owned and j ghost (also stored by proc owning j)
   pair stored once if i,j are both ghost and i < j
------------------------------------------------------------------------- */

void Neighbor::half_nsq_no_newton_ghost_omp(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  const int bitmask = (includegroup) ? group->bitmask[includegroup] : 0;
  const int nall = nlocal + atom->nghost;

  NEIGH_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NEIGH_OMP_SETUP(nall);

  int i,j,n,itype,jtype,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  int **special = atom->special;
  int **nspecial = atom->nspecial;
  int *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int molecular = atom->molecular;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  ipage.reset();

  // loop over owned & ghost atoms, storing neighbors

  for (i = ifrom; i < ito; i++) {

    n = 0;
    neighptr = ipage.vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over remaining atoms, owned and ghost
    // only store pair if i < j
    // stores own/own pairs only once
    // stores own/ghost pairs with owned atom only, on both procs
    // stores ghost/ghost pairs only once
    // no molecular test when i = ghost atom

    if (i < nlocal) {
      for (j = i+1; j < nall; j++) {
        if (includegroup && !(mask[j] & bitmask)) continue;
        jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq[itype][jtype]) {
          if (molecular) {
            which = find_special(special[i],nspecial[i],tag[j]);
            if (which == 0) neighptr[n++] = j;
            else if (domain->minimum_image_check(delx,dely,delz))
              neighptr[n++] = j;
            else if (which > 0) neighptr[n++] = j ^ (which << SBBITS);
          } else neighptr[n++] = j;
        }
      }

    } else {
      for (j = i+1; j < nall; j++) {
        if (includegroup && !(mask[j] & bitmask)) continue;
        jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq[itype][jtype]) neighptr[n++] = j;
      }
    }

    ilist[i] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage.vgot(n);
    if (ipage.status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
  NEIGH_OMP_CLOSE;
  list->inum = atom->nlocal;
  list->gnum = nall - atom->nlocal;
}

/* ----------------------------------------------------------------------
   N^2 / 2 search for neighbor pairs with full Newton's 3rd law
   every pair stored exactly once by some processor
   decision on ghost atoms based on itag,jtag tests
------------------------------------------------------------------------- */

void Neighbor::half_nsq_newton_omp(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  const int bitmask = (includegroup) ? group->bitmask[includegroup] : 0;

  NEIGH_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NEIGH_OMP_SETUP(nlocal);

  int i,j,n,itype,jtype,itag,jtag,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  // loop over each atom, storing neighbors

  int **special = atom->special;
  int **nspecial = atom->nspecial;
  int *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nall = atom->nlocal + atom->nghost;
  int molecular = atom->molecular;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  ipage.reset();

  for (i = ifrom; i < ito; i++) {

    n = 0;
    neighptr = ipage.vget();

    itag = tag[i];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over remaining atoms, owned and ghost
    // itag = jtag is possible for long cutoffs that include images of self

    for (j = i+1; j < nall; j++) {
      if (includegroup && !(mask[j] & bitmask)) continue;

      if (j >= nlocal) {
        jtag = tag[j];
        if (itag > jtag) {
          if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag+jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }
      }

      jtype = type[j];
      if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq[itype][jtype]) {
        if (molecular) {
          which = find_special(special[i],nspecial[i],tag[j]);
          if (which == 0) neighptr[n++] = j;
          else if (domain->minimum_image_check(delx,dely,delz))
            neighptr[n++] = j;
          else if (which > 0) neighptr[n++] = j ^ (which << SBBITS);
        } else neighptr[n++] = j;
      }
    }

    ilist[i] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage.vgot(n);
    if (ipage.status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
  NEIGH_OMP_CLOSE;
  list->inum = nlocal;
}

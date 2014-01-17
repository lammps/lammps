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
#include "domain.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   binned neighbor list construction with partial Newton's 3rd law
   each owned atom i checks own bin and other bins in stencil
   multi-type stencil is itype dependent and is distance checked
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::half_multi_no_newton(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,which,ns;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*s;
  double *cutsq,*distsq;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int molecular = atom->molecular;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms in other bins in stencil including self
    // only store pair if i < j
    // skip if i,j neighbor cutoff is less than bin distance
    // stores own/own pairs only once
    // stores own/ghost pairs on both procs

    ibin = coord2bin(x[i]);
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    cutsq = cutneighsq[itype];
    ns = nstencil_multi[itype];
    for (k = 0; k < ns; k++) {
      for (j = binhead[ibin+s[k]]; j >= 0; j = bins[j]) {
        if (j <= i) continue;
        jtype = type[j];
        if (cutsq[jtype] < distsq[k]) continue;

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
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   multi-type stencil is itype dependent and is distance checked
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::half_multi_newton(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,which,ns;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*s;
  double *cutsq,*distsq;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int molecular = atom->molecular;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over rest of atoms in i's bin, ghosts are at end of linked list
    // if j is owned atom, store it, since j is beyond i in linked list
    // if j is ghost, only store if j coords are "above and to the right" of i

    for (j = bins[i]; j >= 0; j = bins[j]) {
      if (j >= nlocal) {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp) {
          if (x[j][1] < ytmp) continue;
          if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
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

    // loop over all atoms in other bins in stencil, store every pair
    // skip if i,j neighbor cutoff is less than bin distance

    ibin = coord2bin(x[i]);
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    cutsq = cutneighsq[itype];
    ns = nstencil_multi[itype];
    for (k = 0; k < ns; k++) {
      for (j = binhead[ibin+s[k]]; j >= 0; j = bins[j]) {
        jtype = type[j];
        if (cutsq[jtype] < distsq[k]) continue;

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
   binned neighbor list construction with Newton's 3rd law for triclinic
   each owned atom i checks its own bin and other bins in triclinic stencil
   multi-type stencil is itype dependent and is distance checked
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::half_multi_newton_tri(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,which,ns;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*s;
  double *cutsq,*distsq;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int molecular = atom->molecular;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms in bins, including self, in stencil
    // skip if i,j neighbor cutoff is less than bin distance
    // bins below self are excluded from stencil
    // pairs for atoms j "below" i are excluded
    // below = lower z or (equal z and lower y) or (equal zy and lower x)
    //         (equal zyx and j <= i)
    // latter excludes self-self interaction but allows superposed atoms

    ibin = coord2bin(x[i]);
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    cutsq = cutneighsq[itype];
    ns = nstencil_multi[itype];
    for (k = 0; k < ns; k++) {
      for (j = binhead[ibin+s[k]]; j >= 0; j = bins[j]) {
        jtype = type[j];
        if (cutsq[jtype] < distsq[k]) continue;
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp) {
          if (x[j][1] < ytmp) continue;
          if (x[j][1] == ytmp) {
            if (x[j][0] < xtmp) continue;
            if (x[j][0] == xtmp && j <= i) continue;
          }
        }

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

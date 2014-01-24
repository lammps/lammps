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
   N^2 search for all neighbors
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void Neighbor::full_nsq_omp(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  const int bitmask = (includegroup) ? group->bitmask[includegroup] : 0;

  NEIGH_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NEIGH_OMP_SETUP(nlocal);

  int i,j,n,itype,jtype,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int nall = atom->nlocal + atom->nghost;
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

    // loop over all atoms, owned and ghost
    // skip i = j

    for (j = 0; j < nall; j++) {
      if (includegroup && !(mask[j] & bitmask)) continue;
      if (i == j) continue;
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
  list->gnum = 0;
}

/* ----------------------------------------------------------------------
   N^2 search for all neighbors
   include neighbors of ghost atoms, but no "special neighbors" for ghosts
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void Neighbor::full_nsq_ghost_omp(NeighList *list)
{
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  NEIGH_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NEIGH_OMP_SETUP(nall);

  int i,j,n,itype,jtype,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
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

    // loop over all atoms, owned and ghost
    // skip i = j
    // no molecular test when i = ghost atom

    if (i < nlocal) {
      for (j = 0; j < nall; j++) {
        if (i == j) continue;
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
      for (j = 0; j < nall; j++) {
        if (i == j) continue;
        jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighghostsq[itype][jtype]) neighptr[n++] = j;
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
  list->gnum = nall - nlocal;
}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void Neighbor::full_bin_omp(NeighList *list)
{
  // bin owned & ghost atoms

  bin_atoms();

  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;

  NEIGH_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NEIGH_OMP_SETUP(nlocal);

  int i,j,k,n,itype,jtype,ibin,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int molecular = atom->molecular;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;

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

    // loop over all atoms in surrounding bins in stencil including self
    // skip i = j

    ibin = coord2bin(x[i]);

    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
        if (i == j) continue;

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
  list->gnum = 0;
}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   include neighbors of ghost atoms, but no "special neighbors" for ghosts
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void Neighbor::full_bin_ghost_omp(NeighList *list)
{
  // bin owned & ghost atoms

  bin_atoms();

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  NEIGH_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NEIGH_OMP_SETUP(nall);

  int i,j,k,n,itype,jtype,ibin,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int xbin,ybin,zbin,xbin2,ybin2,zbin2;
  int *neighptr;

  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int molecular = atom->molecular;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;
  int **stencilxyz = list->stencilxyz;

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

    // loop over all atoms in surrounding bins in stencil including self
    // when i is a ghost atom, must check if stencil bin is out of bounds
    // skip i = j
    // no molecular test when i = ghost atom

    if (i < nlocal) {
      ibin = coord2bin(x[i]);
      for (k = 0; k < nstencil; k++) {
        for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
          if (i == j) continue;

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
      }

    } else {
      ibin = coord2bin(x[i],xbin,ybin,zbin);
      for (k = 0; k < nstencil; k++) {
        xbin2 = xbin + stencilxyz[k][0];
        ybin2 = ybin + stencilxyz[k][1];
        zbin2 = zbin + stencilxyz[k][2];
        if (xbin2 < 0 || xbin2 >= mbinx ||
            ybin2 < 0 || ybin2 >= mbiny ||
            zbin2 < 0 || zbin2 >= mbinz) continue;
        for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
          if (i == j) continue;

          jtype = type[j];
          if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;

          if (rsq <= cutneighghostsq[itype][jtype]) neighptr[n++] = j;
        }
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
  list->gnum = nall - nlocal;
}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   multi-type stencil is itype dependent and is distance checked
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void Neighbor::full_multi_omp(NeighList *list)
{
  // bin local & ghost atoms

  bin_atoms();

  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;

  NEIGH_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NEIGH_OMP_SETUP(nlocal);

  int i,j,k,n,itype,jtype,ibin,which,ns;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*s;
  double *cutsq,*distsq;

  // loop over each atom, storing neighbors

  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int molecular = atom->molecular;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *nstencil_multi = list->nstencil_multi;
  int **stencil_multi = list->stencil_multi;
  double **distsq_multi = list->distsq_multi;

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  ipage.reset();

  for (i = ifrom; i < ito; i++) {

    n = 0;
    neighptr = ipage.vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms in other bins in stencil, including self
    // skip if i,j neighbor cutoff is less than bin distance
    // skip i = j

    ibin = coord2bin(x[i]);
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    cutsq = cutneighsq[itype];
    ns = nstencil_multi[itype];
    for (k = 0; k < ns; k++) {
      for (j = binhead[ibin+s[k]]; j >= 0; j = bins[j]) {
        jtype = type[j];
        if (cutsq[jtype] < distsq[k]) continue;
        if (i == j) continue;

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

    ilist[i] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage.vgot(n);
    if (ipage.status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
  NEIGH_OMP_CLOSE;
  list->inum = nlocal;
  list->gnum = 0;
}

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

#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   convert atom coords into the ssa active interaction region number
------------------------------------------------------------------------- */

int Neighbor::coord2ssa_airnum(double *x)
{
  int ix, iy, iz;

  ix = iy = iz = 0;
  if (x[2] < domain->sublo[2]) iz = -1;
  if (x[2] > domain->subhi[2]) iz = 1;
  if (x[1] < domain->sublo[1]) iy = -1;
  if (x[1] > domain->subhi[1]) iy = 1;
  if (x[0] < domain->sublo[0]) ix = -1;
  if (x[0] > domain->subhi[0]) ix = 1;

  if(iz < 0) return 0;

  if(iz == 0){
    if( iy<0 ) return 0; // bottom left/middle/right
    if( (iy==0) && (ix<0)  ) return 0; // left atoms
    if( (iy==0) && (ix==0) ) return 1; // Locally owned atoms
    if( (iy==0) && (ix>0)  ) return 3; // Right atoms
    if( (iy>0)  && (ix==0) ) return 2; // Top-middle atoms
    if( (iy>0)  && (ix!=0) ) return 4; // Top-right and top-left atoms
  } else if(iz > 0) {
    if((ix==0) && (iy==0)) return 5; // Back atoms
    if((ix==0) && (iy!=0)) return 6; // Top-back and bottom-back atoms
    if((ix!=0) && (iy==0)) return 7; // Left-back and right-back atoms
    if((ix!=0) && (iy!=0)) return 8; // Back corner atoms
  }

  return 0;
}

/* ----------------------------------------------------------------------
 *    assign owned and ghost atoms their ssa active interaction region numbers
 *    Called in the pre_neighbor and setup_pre_neighbor fix stages
------------------------------------------------------------------------- */

void Neighbor::assign_ssa_airnums()
{
  int i,ibin;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (nall > len_ssa_airnum) {
    len_ssa_airnum = nall;
    memory->destroy(ssa_airnum);
    memory->create(ssa_airnum,len_ssa_airnum,"ssa_airnum");
  }

  // bin in reverse order so linked list will be in forward order

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    for (i = nall-1; i >= nlocal; i--) {
      if (mask[i] & bitmask) {
        ibin = coord2ssa_airnum(x[i]);
      } else {
        ibin = 0;
      }
      ssa_airnum[i] = ibin;
    }
    // All the local excluded atoms are in the zero airnum
    ibin = 0;
    for (i = atom->nlocal-1; i >= atom->nfirst; i--) {
      ssa_airnum[i] = ibin;
    }
  } else {
    for (i = nall-1; i >= nlocal; i--) {
      ibin = coord2ssa_airnum(x[i]);
      ssa_airnum[i] = ibin;
    }
  }

  // All the local included atoms are in the same airnum (#1)
  if (i >= 0) {
    ibin = coord2ssa_airnum(x[i]);
    do {
      ssa_airnum[i] = ibin;
    } while (--i >= 0);
  }
}

/* ----------------------------------------------------------------------
   routines to create a stencil = list of bin offsets
   stencil = bins whose closest corner to central bin is within cutoff
   sx,sy,sz = bin bounds = furthest the stencil could possibly extend
   3d creates xyz stencil, 2d creates xy stencil
   for half list with newton on:
     stencil is bins to the "upper right" of central bin
     stencil does not include self
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_bin_2d_ssa(NeighList *list,
                                   int sx, int sy, int sz)
{
  int i,j;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (j = 0; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (j > 0 || (j == 0 && i > 0))
        if (bin_distance(i,j,0) < cutneighmaxsq)
          stencil[nstencil++] = j*mbinx + i;

  list->nstencil = nstencil;

  // Now include additional bins for AIR ghosts only
  for (j = -sy; j <= 0; j++)
    for (i = -sx; i <= sx; i++) {
      if (j == 0 && i > 0) continue;
      if (bin_distance(i,j,0) < cutneighmaxsq)
        stencil[nstencil++] = j*mbinx + i;
    }

  while (nstencil < list->maxstencil) {
    stencil[nstencil++] = INT_MAX;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_bin_3d_ssa(NeighList *list,
                                   int sx, int sy, int sz)
{
  int i,j,k;
  int *stencil = list->stencil;
  int nstencil = 0;

  for (k = 0; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (k > 0 || j > 0 || (j == 0 && i > 0))
          if (bin_distance(i,j,k) < cutneighmaxsq)
            stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;

  list->nstencil = nstencil;

  // Now include additional bins for AIR ghosts only
  for (k = -sz; k < 0; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq)
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
  k = 0; // skip already included bins at k == 0
  for (j = -sy; j <= 0; j++)
    for (i = -sx; i <= sx; i++) {
      if (j == 0 && i > 0) continue;
      if (bin_distance(i,j,k) < cutneighmaxsq)
        stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
    }

  while (nstencil < list->maxstencil) {
    stencil[nstencil++] = INT_MAX;
  }
}

/* ----------------------------------------------------------------------
   build half list from full list for use by Shardlow Spliting Algorithm
   pair stored once if i,j are both owned and i < j
   if j is ghost, only store if j coords are "above and to the right" of i
   works if full list is a skip list
------------------------------------------------------------------------- */

void Neighbor::half_from_full_newton_ssa(NeighList *list)
{
  int i,j,ii,jj,n,jnum,joriginal;
  int *neighptr,*jlist;

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

    // loop over full neighbor list

    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (j < nlocal) {
        if (i > j) continue;
      } else {
        if (ssa_airnum[j] <= 0) continue;
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
   for Shardlow Spliting Algorithm:
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::half_bin_newton_ssa(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,which,imol,iatom,moltemplate;
  tagint tagprev;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;
  int molecular = atom->molecular;
  if (molecular == 2) moltemplate = 1;
  else moltemplate = 0;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int nstencil = list->nstencil;
  int maxstencil = list->maxstencil;
  int *stencil = list->stencil;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;

/* ----------------------------------------------------------------------
   bin owned and ghost atoms for use by Shardlow Splitting Algorithm
    exclude ghost atoms that are not in the Active Interaction Regions (AIR)
------------------------------------------------------------------------- */

  if (mbins > maxhead_ssa) {
    maxhead_ssa = mbins;
    memory->destroy(gbinhead_ssa);
    memory->destroy(binhead_ssa);
    memory->create(binhead_ssa,maxhead_ssa,"binhead_ssa");
    memory->create(gbinhead_ssa,maxhead_ssa,"gbinhead_ssa");
  }
  for (i = 0; i < mbins; i++) {
    gbinhead_ssa[i] = -1;
    binhead_ssa[i] = -1;
  }

  if (maxbin > maxbin_ssa) {
    maxbin_ssa = maxbin;
    memory->destroy(bins_ssa);
    memory->create(bins_ssa,maxbin_ssa,"bins_ssa");
  }

  // bin in reverse order so linked list will be in forward order

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    for (i = nall-1; i >= nlocal; i--) {
      if (ssa_airnum[i] <= 0) continue; // skip ghost atoms not in AIR
      if (mask[i] & bitmask) {
        ibin = coord2bin(x[i]);
        bins_ssa[i] = gbinhead_ssa[ibin];
        gbinhead_ssa[ibin] = i;
      }
    }
    nlocal = atom->nfirst; // This is important for the code that follows!
  } else {
    for (i = nall-1; i >= nlocal; i--) {
      if (ssa_airnum[i] <= 0) continue; // skip ghost atoms not in AIR
      ibin = coord2bin(x[i]);
      bins_ssa[i] = gbinhead_ssa[ibin];
      gbinhead_ssa[ibin] = i;
    }
  }
  for (i = nlocal-1; i >= 0; i--) {
    ibin = coord2bin(x[i]);
    bins_ssa[i] = binhead_ssa[ibin];
    binhead_ssa[ibin] = i;
  }

  ipage->reset();

  // loop over owned atoms, storing half of the neighbors

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    if (moltemplate) {
      imol = molindex[i];
      iatom = molatom[i];
      tagprev = tag[i] - iatom - 1;
    }

    // loop over rest of local atoms in i's bin
    // just store them, since j is beyond i in linked list

    for (j = bins_ssa[i]; j >= 0; j = bins_ssa[j]) {

      jtype = type[j];
      if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq[itype][jtype]) {
        if (molecular) {
          if (!moltemplate)
            which = find_special(special[i],nspecial[i],tag[j]);
          else if (imol >= 0)
            which = find_special(onemols[imol]->special[iatom],
                                 onemols[imol]->nspecial[iatom],
                                 tag[j]-tagprev);
          else which = 0;
          if (which == 0) neighptr[n++] = j;
          else if (domain->minimum_image_check(delx,dely,delz))
            neighptr[n++] = j;
          else if (which > 0) neighptr[n++] = j ^ (which << SBBITS);
        } else neighptr[n++] = j;
      }
    }

    ibin = coord2bin(x[i]);

    // loop over all local atoms in other bins in "half" stencil
    for (k = 0; k < nstencil; k++) {
      for (j = binhead_ssa[ibin+stencil[k]]; j >= 0; j = bins_ssa[j]) {

        jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq[itype][jtype]) {
          if (molecular) {
            if (!moltemplate)
              which = find_special(special[i],nspecial[i],tag[j]);
            else if (imol >= 0)
              which = find_special(onemols[imol]->special[iatom],
                                   onemols[imol]->nspecial[iatom],
                                   tag[j]-tagprev);
            else which = 0;
            if (which == 0) neighptr[n++] = j;
            else if (domain->minimum_image_check(delx,dely,delz))
              neighptr[n++] = j;
            else if (which > 0) neighptr[n++] = j ^ (which << SBBITS);
          } else neighptr[n++] = j;
        }
      }
    }

    // loop over AIR ghost atoms in all bins in "full" stencil
    // Note: the non-AIR ghost atoms have already been filtered out
    for (k = 0; k < maxstencil; k++) {
      if (stencil[k] > mbins) break; /* Check if ghost stencil bins are exhausted */
      for (j = gbinhead_ssa[ibin+stencil[k]]; j >= 0; j = bins_ssa[j]) {

        jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq[itype][jtype]) {
          if (molecular) {
            if (!moltemplate)
              which = find_special(special[i],nspecial[i],tag[j]);
            else if (imol >= 0)
              which = find_special(onemols[imol]->special[iatom],
                                   onemols[imol]->nspecial[iatom],
                                   tag[j]-tagprev);
            else which = 0;
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

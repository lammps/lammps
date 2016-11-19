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

// space for static variable ssaAIRptr so it
// can be used in qsort's compair function "cmp_ssaAIR()"
static int *ssaAIRptr;

static int cmp_ssaAIR(const void *iptr, const void *jptr)
{
  int i = *((int *) iptr);
  int j = *((int *) jptr);
  if (ssaAIRptr[i] < ssaAIRptr[j]) return -1;
  if (ssaAIRptr[i] > ssaAIRptr[j]) return 1;
  return 0;
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

    // Do a prefix sum on the counts to turn them into indexes.
    list->ndxAIR_ssa[i][0] = AIRct[0];
    for (int ndx = 1; ndx < 8; ++ndx) {
      list->ndxAIR_ssa[i][ndx] = AIRct[ndx] + list->ndxAIR_ssa[i][ndx - 1];
    }
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
  int *ssaAIR = atom->ssaAIR;

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

  if (binatomflag) { /* only false in Neighbor::build_one */
/* ----------------------------------------------------------------------
   bin owned and ghost atoms for use by Shardlow Splitting Algorithm
    exclude ghost atoms that are not in the Active Interaction Regions (AIR)
------------------------------------------------------------------------- */

    if (mbins > list->maxhead_ssa) {
      list->maxhead_ssa = mbins;
      memory->destroy(list->gbinhead_ssa);
      memory->destroy(list->binhead_ssa);
      memory->create(list->binhead_ssa,list->maxhead_ssa,"binhead_ssa");
      memory->create(list->gbinhead_ssa,list->maxhead_ssa,"gbinhead_ssa");
    }
    for (i = 0; i < mbins; i++) {
      list->gbinhead_ssa[i] = -1;
      list->binhead_ssa[i] = -1;
    }

    if (maxbin > list->maxbin_ssa) {
      list->maxbin_ssa = maxbin;
      memory->destroy(list->bins_ssa);
      memory->create(list->bins_ssa,list->maxbin_ssa,"bins_ssa");
    }

    // bin in reverse order so linked list will be in forward order

    if (includegroup) {
      int bitmask = group->bitmask[includegroup];
      for (i = nall-1; i >= nlocal; i--) {
        if (ssaAIR[i] < 2) continue; // skip ghost atoms not in AIR
        if (mask[i] & bitmask) {
          ibin = coord2bin(x[i]);
          list->bins_ssa[i] = list->gbinhead_ssa[ibin];
          list->gbinhead_ssa[ibin] = i;
        }
      }
      nlocal = atom->nfirst; // This is important for the code that follows!
    } else {
      for (i = nall-1; i >= nlocal; i--) {
        if (ssaAIR[i] < 2) continue; // skip ghost atoms not in AIR
        ibin = coord2bin(x[i]);
        list->bins_ssa[i] = list->gbinhead_ssa[ibin];
        list->gbinhead_ssa[ibin] = i;
      }
    }
    for (i = nlocal-1; i >= 0; i--) {
      ibin = coord2bin(x[i]);
      list->bins_ssa[i] = list->binhead_ssa[ibin];
      list->binhead_ssa[ibin] = i;
    }
  } /* else reuse previous binning. See Neighbor::build_one comment. */

  ipage->reset();

  // loop over owned atoms, storing half of the neighbors

  for (i = 0; i < nlocal; i++) {
    int AIRct[8] = { 0 };
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

    for (j = list->bins_ssa[i]; j >= 0; j = list->bins_ssa[j]) {

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
      for (j = list->binhead_ssa[ibin+stencil[k]]; j >= 0; j = list->bins_ssa[j]) {

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
    AIRct[0] = n;

    // loop over AIR ghost atoms in all bins in "full" stencil
    // Note: the non-AIR ghost atoms have already been filtered out
    // That is a significant time savings because of the "full" stencil
    // Note2: only non-pure locals can have ghosts as neighbors
    if (ssaAIR[i] == 1) for (k = 0; k < maxstencil; k++) {
      if (stencil[k] > mbins) break; /* Check if ghost stencil bins are exhausted */
      for (j = list->gbinhead_ssa[ibin+stencil[k]]; j >= 0; j = list->bins_ssa[j]) {

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
            if (which == 0) {
              neighptr[n++] = j;
              ++(AIRct[ssaAIR[j] - 1]);
            } else if (domain->minimum_image_check(delx,dely,delz)) {
              neighptr[n++] = j;
              ++(AIRct[ssaAIR[j] - 1]);
            } else if (which > 0) {
              neighptr[n++] = j ^ (which << SBBITS);
              ++(AIRct[ssaAIR[j] - 1]);
            }
          } else {
            neighptr[n++] = j;
            ++(AIRct[ssaAIR[j] - 1]);
          }
        }
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    // sort the ghosts in the neighbor list by their ssaAIR number
    ssaAIRptr = atom->ssaAIR;
    qsort(&(neighptr[AIRct[0]]), n - AIRct[0], sizeof(int), cmp_ssaAIR);

    // Do a prefix sum on the counts to turn them into indexes.
    list->ndxAIR_ssa[i][0] = AIRct[0];
    for (int ndx = 1; ndx < 8; ++ndx) {
      list->ndxAIR_ssa[i][ndx] = AIRct[ndx] + list->ndxAIR_ssa[i][ndx - 1];
    }
  }

  list->inum = inum;
}

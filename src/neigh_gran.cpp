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
#include "group.h"
#include "fix_shear_history.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   granular particles
   N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
   shear history must be accounted for when a neighbor pair is added
   pair added to list if atoms i and j are both owned and i < j
   pair added if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::granular_nsq_no_newton(NeighList *list)
{
  int i,j,m,n,nn,bitmask;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr,*touchptr;
  double *shearptr;

  NeighList *listgranhistory;
  int *npartner,**partner;
  double ***shearpartner;
  int **firsttouch;
  double **firstshear;
  int **pages_touch;
  double **pages_shear;

  double **x = atom->x;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (includegroup) {
    nlocal = atom->nfirst;
    bitmask = group->bitmask[includegroup];
  }

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;

  FixShearHistory *fix_history = list->fix_history;
  if (fix_history) {
    npartner = fix_history->npartner;
    partner = fix_history->partner;
    shearpartner = fix_history->shearpartner;
    listgranhistory = list->listgranhistory;
    firsttouch = listgranhistory->firstneigh;
    firstshear = listgranhistory->firstdouble;
    pages_touch = listgranhistory->pages;
    pages_shear = listgranhistory->dpages;
  }

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) {
	pages = list->add_pages();
	if (fix_history) {
	  pages_touch = listgranhistory->add_pages();
	  pages_shear = listgranhistory->dpages;
	}
      }
    }

    n = 0;
    neighptr = &pages[npage][npnt];
    if (fix_history) {
      nn = 0;
      touchptr = &pages_touch[npage][npnt];
      shearptr = &pages_shear[npage][3*npnt];
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over remaining atoms, owned and ghost

    for (j = i+1; j < nall; j++) {
      if (includegroup && !(mask[j] & bitmask)) continue;
      if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radius[j];
      cutsq = (radsum+skin) * (radsum+skin);

      if (rsq <= cutsq) {
	neighptr[n] = j;

	if (fix_history) {
	  if (rsq < radsum*radsum) {
	    for (m = 0; m < npartner[i]; m++)
	      if (partner[i][m] == tag[j]) break;
	    if (m < npartner[i]) {
	      touchptr[n] = 1;
	      shearptr[nn++] = shearpartner[i][m][0];
	      shearptr[nn++] = shearpartner[i][m][1];
	      shearptr[nn++] = shearpartner[i][m][2];
	    } else {
	      touchptr[n] = 0;
	      shearptr[nn++] = 0.0;
	      shearptr[nn++] = 0.0;
	      shearptr[nn++] = 0.0;
	    }
	  } else {
	    touchptr[n] = 0;
	    shearptr[nn++] = 0.0;
	    shearptr[nn++] = 0.0;
	    shearptr[nn++] = 0.0;
	  }
	}

	n++;
      }
    }	       

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    if (fix_history) {
      firsttouch[i] = touchptr;
      firstshear[i] = shearptr;
    }
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   granular particles
   N^2 / 2 search for neighbor pairs with full Newton's 3rd law
   no shear history is allowed for this option
   pair added to list if atoms i and j are both owned and i < j
   if j is ghost only me or other proc adds pair
   decision based on itag,jtag tests
------------------------------------------------------------------------- */

void Neighbor::granular_nsq_newton(NeighList *list)
{
  int i,j,n,itag,jtag,bitmask;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr;

  double **x = atom->x;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (includegroup) {
    nlocal = atom->nfirst;
    bitmask = group->bitmask[includegroup];
  }

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }

    n = 0;
    neighptr = &pages[npage][npnt];

    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over remaining atoms, owned and ghost

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

      if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radius[j];
      cutsq = (radsum+skin) * (radsum+skin);
      
      if (rsq <= cutsq) neighptr[n++] = j;
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   granular particles
   binned neighbor list construction with partial Newton's 3rd law
   shear history must be accounted for when a neighbor pair is added
   each owned atom i checks own bin and surrounding bins in non-Newton stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::granular_bin_no_newton(NeighList *list)
{
  int i,j,k,m,n,nn,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr,*touchptr;
  double *shearptr;

  NeighList *listgranhistory;
  int *npartner,**partner;
  double ***shearpartner;
  int **firsttouch;
  double **firstshear;
  int **pages_touch;
  double **pages_shear;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;

  FixShearHistory *fix_history = list->fix_history;
  if (fix_history) {
    npartner = fix_history->npartner;
    partner = fix_history->partner;
    shearpartner = fix_history->shearpartner;
    listgranhistory = list->listgranhistory;
    firsttouch = listgranhistory->firstneigh;
    firstshear = listgranhistory->firstdouble;
    pages_touch = listgranhistory->pages;
    pages_shear = listgranhistory->dpages;
  }

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) {
	pages = list->add_pages();
	if (fix_history) {
	  pages_touch = listgranhistory->add_pages();
	  pages_shear = listgranhistory->dpages;
	}
      }
    }

    n = 0;
    neighptr = &pages[npage][npnt];
    if (fix_history) {
      nn = 0;
      touchptr = &pages_touch[npage][npnt];
      shearptr = &pages_shear[npage][3*npnt];
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    ibin = coord2bin(x[i]);

    // loop over all atoms in surrounding bins in stencil including self
    // only store pair if i < j
    // stores own/own pairs only once
    // stores own/ghost pairs on both procs

    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
	if (j <= i) continue;
	if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	radsum = radi + radius[j];
	cutsq = (radsum+skin) * (radsum+skin);

	if (rsq <= cutsq) {
	  neighptr[n] = j;

	  if (fix_history) {
	    if (rsq < radsum*radsum) {
	      for (m = 0; m < npartner[i]; m++)
		if (partner[i][m] == tag[j]) break;
	      if (m < npartner[i]) {
		touchptr[n] = 1;
		shearptr[nn++] = shearpartner[i][m][0];
		shearptr[nn++] = shearpartner[i][m][1];
		shearptr[nn++] = shearpartner[i][m][2];
	      } else {
		touchptr[n] = 0;
		shearptr[nn++] = 0.0;
		shearptr[nn++] = 0.0;
		shearptr[nn++] = 0.0;
	      }
	    } else {
	      touchptr[n] = 0;
	      shearptr[nn++] = 0.0;
	      shearptr[nn++] = 0.0;
	      shearptr[nn++] = 0.0;
	    }
	  }

	  n++;
	}
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    if (fix_history) {
      firsttouch[i] = touchptr;
      firstshear[i] = shearptr;
    }
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   granular particles
   binned neighbor list construction with full Newton's 3rd law
   no shear history is allowed for this option
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::granular_bin_newton(NeighList *list)
{
  int i,j,k,n,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }

    n = 0;
    neighptr = &pages[npage][npnt];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

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

      if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radius[j];
      cutsq = (radsum+skin) * (radsum+skin);

      if (rsq <= cutsq) neighptr[n++] = j;
    }

    // loop over all atoms in other bins in stencil, store every pair

    ibin = coord2bin(x[i]);
    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
	if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	radsum = radi + radius[j];
	cutsq = (radsum+skin) * (radsum+skin);

	if (rsq <= cutsq) neighptr[n++] = j;
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   granular particles
   binned neighbor list construction with Newton's 3rd law for triclinic
   no shear history is allowed for this option
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::granular_bin_newton_tri(NeighList *list)
{
  int i,j,k,n,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }

    n = 0;
    neighptr = &pages[npage][npnt];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // loop over all atoms in bins in stencil
    // pairs for atoms j "below" i are excluded
    // below = lower z or (equal z and lower y) or (equal zy and lower x)
    //         (equal zyx and j <= i)
    // latter excludes self-self interaction but allows superposed atoms

    ibin = coord2bin(x[i]);
    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
	if (x[j][2] < ztmp) continue;
	if (x[j][2] == ztmp) {
	  if (x[j][1] < ytmp) continue;
	  if (x[j][1] == ytmp) {
	    if (x[j][0] < xtmp) continue;
	    if (x[j][0] == xtmp && j <= i) continue;
	  }
	}

	if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	radsum = radi + radius[j];
	cutsq = (radsum+skin) * (radsum+skin);

	if (rsq <= cutsq) neighptr[n++] = j;
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}

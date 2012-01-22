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
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   multiple respa lists
   N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
   pair added to list if atoms i and j are both owned and i < j
   pair added if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::respa_nsq_no_newton(NeighList *list)
{
  int i,j,n,itype,jtype,n_inner,n_middle,bitmask;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*neighptr_inner,*neighptr_middle;

  // loop over each atom, storing neighbors

  int **special = atom->special;
  int **nspecial = atom->nspecial;
  int *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int molecular = atom->molecular;
  if (includegroup) {
    nlocal = atom->nfirst;
    bitmask = group->bitmask[includegroup];
  }

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;

  NeighList *listinner = list->listinner;
  int *ilist_inner = listinner->ilist;
  int *numneigh_inner = listinner->numneigh;
  int **firstneigh_inner = listinner->firstneigh;
  int **pages_inner = listinner->pages;

  NeighList *listmiddle;
  int *ilist_middle,*numneigh_middle,**firstneigh_middle,**pages_middle;
  int respamiddle = list->respamiddle;
  if (respamiddle) {
    listmiddle = list->listmiddle;
    ilist_middle = listmiddle->ilist;
    numneigh_middle = listmiddle->numneigh;
    firstneigh_middle = listmiddle->firstneigh;
    pages_middle = listmiddle->pages;
  }

  int inum = 0;
  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  int which = 0;

  for (i = 0; i < nlocal; i++) {

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

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over remaining atoms, owned and ghost

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
	  if (which >= 0) neighptr[n++] = j ^ (which << SBBITS);
	} else neighptr[n++] = j;

        if (rsq < cut_inner_sq) {
	  if (which == 0) neighptr_inner[n_inner++] = j;
	  else if (which > 0) neighptr_inner[n_inner++] = j ^ (which << SBBITS);
        }

        if (respamiddle && rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	  if (which == 0) neighptr_middle[n_middle++] = j;
	  else if (which > 0) 
	    neighptr_middle[n_middle++] = j ^ (which << SBBITS);
        }
      }
    }

    ilist[inum] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    ilist_inner[inum] = i;
    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (n_inner > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    if (respamiddle) {
      ilist_middle[inum] = i;
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (n_middle > oneatom)
	error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    }

    inum++;
  }

  list->inum = inum;
  listinner->inum = inum;
  if (respamiddle) listmiddle->inum = inum;
}

/* ----------------------------------------------------------------------
   multiple respa lists
   N^2 / 2 search for neighbor pairs with full Newton's 3rd law
   pair added to list if atoms i and j are both owned and i < j
   if j is ghost only me or other proc adds pair
   decision based on itag,jtag tests
------------------------------------------------------------------------- */

void Neighbor::respa_nsq_newton(NeighList *list)
{
  int i,j,n,itype,jtype,itag,jtag,n_inner,n_middle,bitmask;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*neighptr_inner,*neighptr_middle;

  // loop over each atom, storing neighbors

  int **special = atom->special;
  int **nspecial = atom->nspecial;
  int *tag = atom->tag;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int molecular = atom->molecular;
  if (includegroup) {
    nlocal = atom->nfirst;
    bitmask = group->bitmask[includegroup];
  }

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;

  NeighList *listinner = list->listinner;
  int *ilist_inner = listinner->ilist;
  int *numneigh_inner = listinner->numneigh;
  int **firstneigh_inner = listinner->firstneigh;
  int **pages_inner = listinner->pages;

  NeighList *listmiddle;
  int *ilist_middle,*numneigh_middle,**firstneigh_middle,**pages_middle;
  int respamiddle = list->respamiddle;
  if (respamiddle) {
    listmiddle = list->listmiddle;
    ilist_middle = listmiddle->ilist;
    numneigh_middle = listmiddle->numneigh;
    firstneigh_middle = listmiddle->firstneigh;
    pages_middle = listmiddle->pages;
  }

  int inum = 0;
  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  int which = 0;

  for (i = 0; i < nlocal; i++) {

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

    itag = tag[i];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

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

      jtype = type[j];
      if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq[itype][jtype]) {
	if (molecular) {
	  which = find_special(special[i],nspecial[i],tag[j]);
	  if (which >= 0) neighptr[n++] = j ^ (which << SBBITS);
	} else neighptr[n++] = j;

        if (rsq < cut_inner_sq) {
	  if (which == 0) neighptr_inner[n_inner++] = j;
	  else if (which > 0) neighptr_inner[n_inner++] = j ^ (which << SBBITS);
        }

        if (respamiddle && 
	    rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	  if (which == 0) neighptr_middle[n_middle++] = j;
	  else if (which > 0) 
	    neighptr_middle[n_middle++] = j ^ (which << SBBITS);
        }
      }
    }

    ilist[inum] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    ilist_inner[inum] = i;
    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (n_inner > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    if (respamiddle) {
      ilist_middle[inum] = i;
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (n_middle > oneatom)
	error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    }

    inum++;
  }

  list->inum = inum;
  listinner->inum = inum;
  if (respamiddle) listmiddle->inum = inum;
}

/* ----------------------------------------------------------------------
   multiple respa lists
   binned neighbor list construction with partial Newton's 3rd law
   each owned atom i checks own bin and surrounding bins in non-Newton stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::respa_bin_no_newton(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,n_inner,n_middle;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*neighptr_inner,*neighptr_middle;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  int **special = atom->special;
  int **nspecial = atom->nspecial;
  int *tag = atom->tag;

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
  int **pages = list->pages;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;

  NeighList *listinner = list->listinner;
  int *ilist_inner = listinner->ilist;
  int *numneigh_inner = listinner->numneigh;
  int **firstneigh_inner = listinner->firstneigh;
  int **pages_inner = listinner->pages;

  NeighList *listmiddle;
  int *ilist_middle,*numneigh_middle,**firstneigh_middle,**pages_middle;
  int respamiddle = list->respamiddle;
  if (respamiddle) {
    listmiddle = list->listmiddle;
    ilist_middle = listmiddle->ilist;
    numneigh_middle = listmiddle->numneigh;
    firstneigh_middle = listmiddle->firstneigh;
    pages_middle = listmiddle->pages;
  }

  int inum = 0;
  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  int which = 0;

  for (i = 0; i < nlocal; i++) {

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

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    ibin = coord2bin(x[i]);

    // loop over all atoms in surrounding bins in stencil including self
    // only store pair if i < j
    // stores own/own pairs only once
    // stores own/ghost pairs on both procs

    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
	if (j <= i) continue;

	jtype = type[j];
	if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq <= cutneighsq[itype][jtype]) {
	  if (molecular) {
	    which = find_special(special[i],nspecial[i],tag[j]);
	    if (which >= 0) neighptr[n++] = j ^ (which << SBBITS);
	  } else neighptr[n++] = j;

	  if (rsq < cut_inner_sq) {
	    if (which == 0) neighptr_inner[n_inner++] = j;
	    else if (which > 0) 
	      neighptr_inner[n_inner++] = j ^ (which << SBBITS);
	  }

	  if (respamiddle && 
	      rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	    if (which == 0) neighptr_middle[n_middle++] = j;
	    else if (which > 0) 
	      neighptr_middle[n_middle++] = j ^ (which << SBBITS);
	  }
	}
      }
    }

    ilist[inum] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    ilist_inner[inum] = i;
    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (n_inner > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    if (respamiddle) {
      ilist_middle[inum] = i;
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (n_middle > oneatom)
	error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    }

    inum++;
  }

  list->inum = inum;
  listinner->inum = inum;
  if (respamiddle) listmiddle->inum = inum;
}
      
/* ----------------------------------------------------------------------
   multiple respa lists
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::respa_bin_newton(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,n_inner,n_middle;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*neighptr_inner,*neighptr_middle;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  int **special = atom->special;
  int **nspecial = atom->nspecial;
  int *tag = atom->tag;

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
  int **pages = list->pages;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;

  NeighList *listinner = list->listinner;
  int *ilist_inner = listinner->ilist;
  int *numneigh_inner = listinner->numneigh;
  int **firstneigh_inner = listinner->firstneigh;
  int **pages_inner = listinner->pages;

  NeighList *listmiddle;
  int *ilist_middle,*numneigh_middle,**firstneigh_middle,**pages_middle;
  int respamiddle = list->respamiddle;
  if (respamiddle) {
    listmiddle = list->listmiddle;
    ilist_middle = listmiddle->ilist;
    numneigh_middle = listmiddle->numneigh;
    firstneigh_middle = listmiddle->firstneigh;
    pages_middle = listmiddle->pages;
  }

  int inum = 0;
  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  int which = 0;

  for (i = 0; i < nlocal; i++) {

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
	  if (which >= 0) neighptr[n++] = j ^ (which << SBBITS);
	} else neighptr[n++] = j;

        if (rsq < cut_inner_sq) {
	  if (which == 0) neighptr_inner[n_inner++] = j;
	  else if (which > 0) neighptr_inner[n_inner++] = j ^ (which << SBBITS);
        }

        if (respamiddle && 
	    rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	  if (which == 0) neighptr_middle[n_middle++] = j;
	  else if (which > 0) 
	    neighptr_middle[n_middle++] = j ^ (which << SBBITS);
        }
      }
    }

    // loop over all atoms in other bins in stencil, store every pair

    ibin = coord2bin(x[i]);
    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
	jtype = type[j];
	if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq <= cutneighsq[itype][jtype]) {
	  if (molecular) {
	    which = find_special(special[i],nspecial[i],tag[j]);
	    if (which >= 0) neighptr[n++] = j ^ (which << SBBITS);
	  } else neighptr[n++] = j;

	  if (rsq < cut_inner_sq) {
	    if (which == 0) neighptr_inner[n_inner++] = j;
	    else if (which > 0) 
	      neighptr_inner[n_inner++] = j ^ (which << SBBITS);
	  }

	  if (respamiddle && 
	      rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	    if (which == 0) neighptr_middle[n_middle++] = j;
	    else if (which > 0) 
	      neighptr_middle[n_middle++] = j ^ (which << SBBITS);
	  }
	}
      }
    }

    ilist[inum] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    ilist_inner[inum] = i;
    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (n_inner > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    if (respamiddle) {
      ilist_middle[inum] = i;
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (n_middle > oneatom)
	error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    }

    inum++;
  }

  list->inum = inum;
  listinner->inum = inum;
  if (respamiddle) listmiddle->inum = inum;
}

/* ----------------------------------------------------------------------
   multiple respa lists
   binned neighbor list construction with Newton's 3rd law for triclinic
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::respa_bin_newton_tri(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,n_inner,n_middle;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*neighptr_inner,*neighptr_middle;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  int **special = atom->special;
  int **nspecial = atom->nspecial;
  int *tag = atom->tag;

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
  int **pages = list->pages;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;

  NeighList *listinner = list->listinner;
  int *ilist_inner = listinner->ilist;
  int *numneigh_inner = listinner->numneigh;
  int **firstneigh_inner = listinner->firstneigh;
  int **pages_inner = listinner->pages;

  NeighList *listmiddle;
  int *ilist_middle,*numneigh_middle,**firstneigh_middle,**pages_middle;
  int respamiddle = list->respamiddle;
  if (respamiddle) {
    listmiddle = list->listmiddle;
    ilist_middle = listmiddle->ilist;
    numneigh_middle = listmiddle->numneigh;
    firstneigh_middle = listmiddle->firstneigh;
    pages_middle = listmiddle->pages;
  }

  int inum = 0;
  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  int which = 0;

  for (i = 0; i < nlocal; i++) {

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

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

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

	jtype = type[j];
	if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq <= cutneighsq[itype][jtype]) {
	  if (molecular) {
	    which = find_special(special[i],nspecial[i],tag[j]);
	    if (which >= 0) neighptr[n++] = j ^ (which << SBBITS);
	  } else neighptr[n++] = j;

	  if (rsq < cut_inner_sq) {
	    if (which == 0) neighptr_inner[n_inner++] = j;
	    else if (which > 0) 
	      neighptr_inner[n_inner++] = j ^ (which << SBBITS);
	  }

	  if (respamiddle &&
	      rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	    if (which == 0) neighptr_middle[n_middle++] = j;
	    else if (which > 0) 
	      neighptr_middle[n_middle++] = j ^ (which << SBBITS);
	  }
	}
      }
    }

    ilist[inum] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    ilist_inner[inum] = i;
    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (n_inner > oneatom)
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    if (respamiddle) {
      ilist_middle[inum] = i;
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (n_middle > oneatom)
	error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    }

    inum++;
  }

  list->inum = inum;
  listinner->inum = inum;
  if (respamiddle) listmiddle->inum = inum;
}

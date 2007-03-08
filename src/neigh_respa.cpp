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
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   multiple respa lists
   N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
   pair added to list if atoms i and j are both owned and i < j
   pair added if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::respa_nsq_no_newton()
{
  int i,j,itype,jtype,which;
  int n_inner,n_middle,n;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr_inner;
  int *neighptr_middle;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int molecular = atom->molecular;

  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) add_pages(npage);
    }
    neighptr = &pages[npage][npnt];
    n = 0;

    if (pgsize - npnt_inner < oneatom) {
      npnt_inner = 0;
      npage_inner++;
      if (npage_inner == maxpage_inner) add_pages_inner(npage_inner);
    }
    neighptr_inner = &pages_inner[npage_inner][npnt_inner];
    n_inner = 0;

    if (respa == 2) {
      if (pgsize - npnt_middle < oneatom) {
	npnt_middle = 0;
	npage_middle++;
	if (npage_middle == maxpage_middle) add_pages_middle(npage_middle);
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
      if (exclude && exclusion(i,j,type,mask,molecule)) continue;

      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq[itype][jtype]) {
	if (molecular) which = find_special(i,j);
	else which = 0;
	if (which == 0) neighptr[n++] = j;
	else if (which > 0) neighptr[n++] = which*nall + j;

        if (rsq < cut_inner_sq) {
	  if (which == 0) neighptr_inner[n_inner++] = j;
	  else if (which > 0) neighptr_inner[n_inner++] = which*nall + j;
        }

        if (respa == 2 && rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	  if (which == 0) neighptr_middle[n_middle++] = j;
	  else if (which > 0) neighptr_middle[n_middle++] = which*nall + j;
        }
      }
    }

    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (npnt_inner >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    if (respa == 2) {
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (npnt_middle >= pgsize)
	error->one("Neighbor list overflow, boost neigh_modify one or page");
    }
  }
}

/* ----------------------------------------------------------------------
   multiple respa lists
   N^2 / 2 search for neighbor pairs with full Newton's 3rd law
   pair added to list if atoms i and j are both owned and i < j
   if j is ghost only me or other proc adds pair
   decision based on itag,jtag tests
------------------------------------------------------------------------- */

void Neighbor::respa_nsq_newton()
{
  int i,j,itype,jtype,itag,jtag,which;
  int n_inner,n_middle,n;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr_inner;
  int *neighptr_middle;
  int *neighptr;

  double **x = atom->x;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int molecular = atom->molecular;

  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) add_pages(npage);
    }
    neighptr = &pages[npage][npnt];
    n = 0;

    if (pgsize - npnt_inner < oneatom) {
      npnt_inner = 0;
      npage_inner++;
      if (npage_inner == maxpage_inner) add_pages_inner(npage_inner);
    }
    neighptr_inner = &pages_inner[npage_inner][npnt_inner];
    n_inner = 0;

    if (respa == 2) {
      if (pgsize - npnt_middle < oneatom) {
	npnt_middle = 0;
	npage_middle++;
	if (npage_middle == maxpage_middle) add_pages_middle(npage_middle);
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
      if (j >= nlocal) {
	jtag = tag[j];
	if (itag > jtag) {
	  if ((itag+jtag) % 2 == 0) continue;
	} else if (itag < jtag) {
	  if ((itag+jtag) % 2 == 1) continue;
	} else {
	  if (x[j][2] < ztmp) continue;
	  else if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	  else if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp)
	    continue;
	}
      }

      if (exclude && exclusion(i,j,type,mask,molecule)) continue;

      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq[itype][jtype]) {
	if (molecular) which = find_special(i,j);
	else which = 0;
	if (which == 0) neighptr[n++] = j;
	else if (which > 0) neighptr[n++] = which*nall + j;

        if (rsq < cut_inner_sq) {
	  if (which == 0) neighptr_inner[n_inner++] = j;
	  else if (which > 0) neighptr_inner[n_inner++] = which*nall + j;
        }

        if (respa == 2 && rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	  if (which == 0) neighptr_middle[n_middle++] = j;
	  else if (which > 0) neighptr_middle[n_middle++] = which*nall + j;
        }
      }
    }

    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (npnt_inner >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    if (respa == 2) {
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (npnt_middle >= pgsize)
	error->one("Neighbor list overflow, boost neigh_modify one or page");
    }
  }
}

/* ----------------------------------------------------------------------
   multiple respa lists
   binned neighbor list construction with partial Newton's 3rd law
   each owned atom i checks own bin and surrounding bins in non-Newton stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::respa_bin_no_newton()
{
  int i,j,k,itype,jtype,ibin,which;
  int n_inner,n_middle,n;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr_inner;
  int *neighptr_middle;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int molecular = atom->molecular;

  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) add_pages(npage);
    }
    neighptr = &pages[npage][npnt];
    n = 0;

    if (pgsize - npnt_inner < oneatom) {
      npnt_inner = 0;
      npage_inner++;
      if (npage_inner == maxpage_inner) add_pages_inner(npage_inner);
    }
    neighptr_inner = &pages_inner[npage_inner][npnt_inner];
    n_inner = 0;

    if (respa == 2) {
      if (pgsize - npnt_middle < oneatom) {
	npnt_middle = 0;
	npage_middle++;
	if (npage_middle == maxpage_middle) add_pages_middle(npage_middle);
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
	if (exclude && exclusion(i,j,type,mask,molecule)) continue;

	jtype = type[j];
	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq <= cutneighsq[itype][jtype]) {
	  if (molecular) which = find_special(i,j);
	  else which = 0;
	  if (which == 0) neighptr[n++] = j;
	  else if (which > 0) neighptr[n++] = which*nall + j;

	  if (rsq < cut_inner_sq) {
	    if (which == 0) neighptr_inner[n_inner++] = j;
	    else if (which > 0) neighptr_inner[n_inner++] = which*nall + j;
	  }

	  if (respa == 2 && rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	    if (which == 0) neighptr_middle[n_middle++] = j;
	    else if (which > 0) neighptr_middle[n_middle++] = which*nall + j;
	  }
	}
      }
    }

    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (npnt_inner >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    if (respa == 2) {
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (npnt_middle >= pgsize)
	error->one("Neighbor list overflow, boost neigh_modify one or page");
    }
  }
}
      
/* ----------------------------------------------------------------------
   multiple respa lists
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::respa_bin_newton()
{
  int i,j,k,itype,jtype,ibin,which;
  int n_inner,n_middle,n;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr_inner;
  int *neighptr_middle;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int molecular = atom->molecular;

  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) add_pages(npage);
    }
    neighptr = &pages[npage][npnt];
    n = 0;

    if (pgsize - npnt_inner < oneatom) {
      npnt_inner = 0;
      npage_inner++;
      if (npage_inner == maxpage_inner) add_pages_inner(npage_inner);
    }
    neighptr_inner = &pages_inner[npage_inner][npnt_inner];
    n_inner = 0;

    if (respa == 2) {
      if (pgsize - npnt_middle < oneatom) {
	npnt_middle = 0;
	npage_middle++;
	if (npage_middle == maxpage_middle) add_pages_middle(npage_middle);
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
	if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
	if (exclude && exclusion(i,j,type,mask,molecule)) continue;
      }

      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq[itype][jtype]) {
	if (molecular) which = find_special(i,j);
	else which = 0;
	if (which == 0) neighptr[n++] = j;
	else if (which > 0) neighptr[n++] = which*nall + j;

        if (rsq < cut_inner_sq) {
	  if (which == 0) neighptr_inner[n_inner++] = j;
	  else if (which > 0) neighptr_inner[n_inner++] = which*nall + j;
        }

        if (respa == 2 && rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	  if (which == 0) neighptr_middle[n_middle++] = j;
	  else if (which > 0) neighptr_middle[n_middle++] = which*nall + j;
        }
      }
    }

    // loop over all atoms in other bins in stencil, store every pair

    ibin = coord2bin(x[i]);
    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
	if (exclude && exclusion(i,j,type,mask,molecule)) continue;

	jtype = type[j];
	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq <= cutneighsq[itype][jtype]) {
	  if (molecular) which = find_special(i,j);
	  else which = 0;
	  if (which == 0) neighptr[n++] = j;
	  else if (which > 0) neighptr[n++] = which*nall + j;

	  if (rsq < cut_inner_sq) {
	    if (which == 0) neighptr_inner[n_inner++] = j;
	    else if (which > 0) neighptr_inner[n_inner++] = which*nall + j;
	  }

	  if (respa == 2 && rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	    if (which == 0) neighptr_middle[n_middle++] = j;
	    else if (which > 0) neighptr_middle[n_middle++] = which*nall + j;
	  }
	}
      }
    }

    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (npnt_inner >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    if (respa == 2) {
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (npnt_middle >= pgsize)
	error->one("Neighbor list overflow, boost neigh_modify one or page");
    }
  }
}

/* ----------------------------------------------------------------------
   multiple respa lists
   binned neighbor list construction with Newton's 3rd law for triclinic
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::respa_bin_newton_tri()
{
  int i,j,k,itype,jtype,ibin,which;
  int n_inner,n_middle,n;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr_inner;
  int *neighptr_middle;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int molecular = atom->molecular;

  int npage = 0;
  int npnt = 0;
  int npage_inner = 0;
  int npnt_inner = 0;
  int npage_middle = 0;
  int npnt_middle = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) add_pages(npage);
    }
    neighptr = &pages[npage][npnt];
    n = 0;

    if (pgsize - npnt_inner < oneatom) {
      npnt_inner = 0;
      npage_inner++;
      if (npage_inner == maxpage_inner) add_pages_inner(npage_inner);
    }
    neighptr_inner = &pages_inner[npage_inner][npnt_inner];
    n_inner = 0;

    if (respa == 2) {
      if (pgsize - npnt_middle < oneatom) {
	npnt_middle = 0;
	npage_middle++;
	if (npage_middle == maxpage_middle) add_pages_middle(npage_middle);
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
    // below = lower z or (equal z and lower y) or (equal zy and <= x)
    // this excludes self-self interaction

    ibin = coord2bin(x[i]);
    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
	if (x[j][2] < ztmp) continue;
	if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] <= xtmp) continue;
	if (exclude && exclusion(i,j,type,mask,molecule)) continue;

	jtype = type[j];
	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq <= cutneighsq[itype][jtype]) {
	  if (molecular) which = find_special(i,j);
	  else which = 0;
	  if (which == 0) neighptr[n++] = j;
	  else if (which > 0) neighptr[n++] = which*nall + j;

	  if (rsq < cut_inner_sq) {
	    if (which == 0) neighptr_inner[n_inner++] = j;
	    else if (which > 0) neighptr_inner[n_inner++] = which*nall + j;
	  }

	  if (respa == 2 && rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
	    if (which == 0) neighptr_middle[n_middle++] = j;
	    else if (which > 0) neighptr_middle[n_middle++] = which*nall + j;
	  }
	}
      }
    }

    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    npnt_inner += n_inner;
    if (npnt_inner >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");

    if (respa == 2) {
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      npnt_middle += n_middle;
      if (npnt_middle >= pgsize)
	error->one("Neighbor list overflow, boost neigh_modify one or page");
    }
  }
}

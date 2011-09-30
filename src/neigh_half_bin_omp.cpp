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
#include "comm.h"
#include "error.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::half_bin_newton_omp(NeighList *list) 
{

  // bin local & ghost atoms

  bin_atoms();

  const int nthreads = comm->nthreads;
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;

  // make sure we have at least one page for each thread
  if (nthreads > list->maxpage) {
    list->add_pages(nthreads - list->maxpage);
  }
  int inum = 0;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list,inum)
#endif
  {

    // get thread id and then assign each thread a fixed chunk of atoms
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    const int idelta = 1 + nlocal/nthreads;
    const int ifrom = tid*idelta;
    int ito   = ifrom + idelta;
    if (ito > nlocal)
      ito = nlocal;

    int i,j,k,n,itype,jtype,ibin,which;
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
    int molecular = atom->molecular;

    int *ilist = list->ilist;
    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;
    int nstencil = list->nstencil;
    int *stencil = list->stencil;
  
    // each thread works on its own page
    int npage = tid;
    int npnt = 0;

    for (i = ifrom; i < ito; i++) {

      if (pgsize - npnt < oneatom) {
	npnt = 0;
	npage += nthreads;

	// only one thread at a time may check
	// whether we need new neighbor list pages
	// and then add to them.
#if defined(_OPENMP)
#pragma omp critical
#endif
	if (npage >= list->maxpage) {
	  list->add_pages(nthreads);
	}
      }

      neighptr = &(list->pages[npage][npnt]);
      n = 0;

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
	  }
	}
      }

#if defined(_OPENMP)
#pragma omp critical
#endif
      ilist[inum++] = i;

      firstneigh[i] = neighptr;
      numneigh[i] = n;
      npnt += n;
      if (n > oneatom || npnt >= pgsize) {
	error->one(FLERR,"Neighbor list overflow, boost neigh_modify one or page");
      }
    }
    list->inum = inum;
  }
}


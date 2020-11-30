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

#include "omp_compat.h"
#include "npair_half_size_multi2_newton_tri_omp.h"
#include "npair_omp.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfSizeMulti2NewtonTriOmp::NPairHalfSizeMulti2NewtonTriOmp(LAMMPS *lmp) :
  NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction with Newton's 3rd law for triclinic
   each owned atom i checks its own bin and other bins in triclinic stencil
   multi-type stencil is itype dependent and is distance checked
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void NPairHalfSizeMulti2NewtonTriOmp::build(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;  
  const int history = list->history;
  const int mask_history = 3 << SBBITS;

  NPAIR_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(list)
#endif
  NPAIR_OMP_SETUP(nlocal);

  int i,j,k,n,itype,jtype,ibin,jbin,ns;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutdistsq;
  int *neighptr,*s;
  int js;

  // loop over each atom, storing neighbors

  double **x = atom->x;
  double *radius = atom->radius;  
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

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
    radi = radius[i];

    // own type: loop over atoms ahead in bin, including ghosts at end of list
    // if j is owned atom, store by virtue of being ahead of i in list
    // if j is ghost, store if x[j] "above and to right of" x[i]

    ibin = atom2bin_multi2[itype][i];

    for (jtype = 1; jtype <= atom->ntypes; jtype++) {

      if (itype == jtype) {

	    // loop over all atoms in other bins in stencil, store every pair
	    // skip if i,j neighbor cutoff is less than bin distance
        
	    s = stencil_multi2[itype][itype];
	    ns = nstencil_multi2[itype][itype];
	    for (k = 0; k < ns; k++) {
	      js = binhead_multi2[itype][ibin + s[k]];
	      for (j = js; j >= 0; j = bins_multi2[itype][j]) {
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
	        radsum = radi + radius[j];
	        cutdistsq = (radsum+skin) * (radsum+skin);
        
	        if (rsq <= cutdistsq) {
	          if (history && rsq < radsum*radsum) 
	    	    neighptr[n++] = j ^ mask_history;
	          else
	    	    neighptr[n++] = j;
	        }
	      }
	    }
      } else {
        // smaller -> larger: locate i in the jtype bin structure
        
	    jbin = coord2bin(x[i], jtype);
        
	    s = stencil_multi2[itype][jtype];
	    ns = nstencil_multi2[itype][jtype];
	    for (k = 0; k < ns; k++) {
	      js = binhead_multi2[jtype][jbin + s[k]];
	      for (j = js; j >= 0; j = bins_multi2[jtype][j]) {
                    
            // if same size, use half stencil            
            if(cutneighsq[itype][itype] == cutneighsq[jtype][jtype]){
              if (x[j][2] < ztmp) continue;
              if (x[j][2] == ztmp) {
                if (x[j][1] < ytmp) continue;
                if (x[j][1] == ytmp) {
                  if (x[j][0] < xtmp) continue;
                  if (x[j][0] == xtmp && j <= i) continue;
                }
              }                
            }
        
	        if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;
        
	        delx = xtmp - x[j][0];
	        dely = ytmp - x[j][1];
	        delz = ztmp - x[j][2];
	        rsq = delx*delx + dely*dely + delz*delz;
	        radsum = radi + radius[j];
	        cutdistsq = (radsum+skin) * (radsum+skin);
        
	        if (rsq <= cutdistsq) {
	          if (history && rsq < radsum*radsum) 
	    	    neighptr[n++] = j ^ mask_history;
	          else
	    	    neighptr[n++] = j;
	        }
	      }
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
  NPAIR_OMP_CLOSE;
  list->inum = nlocal;
}

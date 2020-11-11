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

#include "nstencil_half_multi2_3d_tri.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "nbin.h"
#include "memory.h"
#include "atom.h"
#include <math.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfMulti23dTri::NStencilHalfMulti23d(LAMMPS *lmp) :
  NStencil(lmp) {}

/* ---------------------------------------------------------------------- */

void NStencilHalfMulti23d::set_stencil_properties()
{
  int n = atom->ntypes;
  int i, j;
  
  // like -> like => use half stencil
  for (i = 1; i <= n; i++) {
    stencil_half[i][i] = 1;
    stencil_skip[i][i] = 0;
    stencil_bin_type[i][i] = i;
    stencil_cut[i][i] = sqrt(cutneighsq[i][i]);
  }

  // Cross types: use full stencil, looking one way through hierarchy
  // smaller -> larger => use full stencil in larger bin
  // larger -> smaller => no nstecil required
  // If cut offs are same, use existing type-type stencil

  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      if(i == j) continue;
      if(cuttypesq[i] > cuttypesq[j]) continue;

      stencil_skip[i][j] = 0;
      stencil_cut[i][j] = sqrt(cutneighsq[i][j]);          
      
      if(cuttypesq[i] == cuttypesq[j]){
        stencil_half[i][j] = 1;
        stencil_bin_type[i][j] = i;
      } else {
        stencil_half[i][j] = 0;
        stencil_bin_type[i][j] = j;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create stencils based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilHalfMulti23d::create()
{
  int itype, jtype, i, j, k, ns;
  int n = atom->ntypes;
  double cutsq;
  
  
  for (itype = 1; itype <= n; itype++) {
    for (jtype = 1; jtype <= n; jtype++) {
      if (stencil_skip[itype][jtype]) continue;
      
      ns = 0;
      
      sx = sx_multi2[itype][jtype];
      sy = sy_multi2[itype][jtype];
      sz = sz_multi2[itype][jtype];
      
      mbinx = mbinx_multi2[itype][jtype];
      mbiny = mbiny_multi2[itype][jtype];
      mbinz = mbinz_multi2[itype][jtype];
      
      cutsq = stencil_cut[itype][jtype];
      
      if (stencil_half[itype][jtype]) {
        for (k = 0; k <= sz; k++)
          for (j = -sy; j <= sy; j++)
            for (i = -sx; i <= sx; i++)
              if (bin_distance(i,j,k) < cutsq)
	            stencil_type[itype][jtype][ns++] = 
                        k*mbiny*mbinx + j*mbinx + i;
      } else {
        for (k = -sz; k <= sz; k++)
          for (j = -sy; j <= sy; j++)
            for (i = -sx; i <= sx; i++)
	          if (bin_distance(i,j,k) < cutsq)
	            stencil_type[itype][jtype][ns++] = 
                        k*mbiny*mbinx + j*mbinx + i;
      }
      
      nstencil_multi2[itype][jtype] = ns;
    }
  }
}

/* ---------------------------------------------------------------------- */

// KS To superclass
void NStencilHalfMulti23d::copy_bin_info_bytype(int itype) {

  mbinx = nb->mbinx_type[itype];
  mbiny = nb->mbiny_type[itype];
  mbinz = nb->mbinz_type[itype];
  binsizex = nb->binsizex_type[itype];
  binsizey = nb->binsizey_type[itype];
  binsizez = nb->binsizez_type[itype];
  bininvx = nb->bininvx_type[itype];
  bininvy = nb->bininvy_type[itype];
  bininvz = nb->bininvz_type[itype];
}

/* ---------------------------------------------------------------------- */

// KS To superclass?
int NStencilHalfMulti23d::copy_neigh_info_bytype(int itype) {

  cutneighmaxsq = neighbor->cutneighsq[itype][itype];
  cutneighmax   = sqrt(cutneighmaxsq);
  cuttypesq     = neighbor->cuttypesq;

  // sx,sy,sz = max range of stencil in each dim
  // smax = max possible size of entire 3d stencil
  // stencil will be empty if cutneighmax = 0.0

  sx = static_cast<int> (cutneighmax*bininvx);
  if (sx*binsizex < cutneighmax) sx++;
  sy = static_cast<int> (cutneighmax*bininvy);
  if (sy*binsizey < cutneighmax) sy++;
  sz = static_cast<int> (cutneighmax*bininvz);
  if (sz*binsizez < cutneighmax) sz++;

  return ((2*sx+1) * (2*sy+1) * (2*sz+1));
}

/* ---------------------------------------------------------------------- */

void NStencilHalfMulti23d::create_setup()
{

  int itype, jtype;
  int maxtypes;
  int smax;

  // maxstencil_type to superclass?
  maxtypes = atom->ntypes;

  if (maxstencil_type == NULL) {
    memory->create(maxstencil_type, maxtypes+1, maxtypes+1, "maxstencil_type");
    memory->create(nstencil_type, maxtypes+1, maxtypes+1, "nstencil_type");
    stencil_type = new int**[maxtypes+1]();
    for (itype = 1; itype <= maxtypes; ++itype) {
      stencil_type[itype] = new int*[maxtypes+1]();
      for (jtype = 1; jtype <= maxtypes; ++jtype) {
	maxstencil_type[itype][jtype] = 0;
	nstencil_type[itype][jtype] = 0;
      }
    }
  }

  // like -> like => use standard Newton stencil at bin 

  for (itype = 1; itype <= maxtypes; ++itype) {
    copy_bin_info_bytype(itype);
    smax = copy_neigh_info_bytype(itype);
    if (smax > maxstencil_type[itype][itype]) {
      maxstencil_type[itype][itype] = smax;
      memory->destroy(stencil_type[itype][itype]);
      memory->create(stencil_type[itype][itype], smax,
		     "NStencilHalfBytypeNewton::create_steup() stencil");
    }
    create_newton(itype, itype, cutneighmaxsq);
  }

  // Cross types: "Newton on" reached by using Newton off stencil and
  // looking one way through hierarchy
  // smaller -> larger => use Newton off stencil in larger bin
  // larger -> smaller => no nstecil required
  // If cut offs are same, use existing type-type stencil

  for (itype = 1; itype <= maxtypes; ++itype) {
    for (jtype = 1; jtype <= maxtypes; ++jtype) {
      if (itype == jtype) continue;
      if (cuttypesq[itype] == cuttypesq[jtype]) {
	nstencil_type[itype][jtype] = nstencil_type[jtype][jtype];
	stencil_type[itype][jtype]  = stencil_type[jtype][jtype];
      }
      else if (cuttypesq[itype] < cuttypesq[jtype]) {
	copy_bin_info_bytype(jtype);

	cutneighmaxsq = cuttypesq[jtype];
	cutneighmax   = sqrt(cutneighmaxsq);
	sx = static_cast<int> (cutneighmax*bininvx);
	if (sx*binsizex < cutneighmax) sx++;
	sy = static_cast<int> (cutneighmax*bininvy);
	if (sy*binsizey < cutneighmax) sy++;
	sz = static_cast<int> (cutneighmax*bininvz);
	if (sz*binsizez < cutneighmax) sz++;

	smax = (2*sx+1) * (2*sy+1) * (2*sz+1);
	if (smax > maxstencil_type[itype][jtype]) {
	  maxstencil_type[itype][jtype] = smax;
	  memory->destroy(stencil_type[itype][jtype]);
	  memory->create(stencil_type[itype][jtype], smax, "stencil_type[]");
	}
	create_newtoff(itype, jtype, cuttypesq[jtype]);
      }
    }
  }
}

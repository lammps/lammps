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

#include "nstencil_half_bytype_3d_newtoff.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "nbin.h"
#include "memory.h"
#include "atom.h"

#include <math.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfBytype3dNewtoff::NStencilHalfBytype3dNewtoff(LAMMPS *lmp) :
  NStencil(lmp)
{
  maxstencil_type = NULL;
}

NStencilHalfBytype3dNewtoff::~NStencilHalfBytype3dNewtoff() {

  memory->destroy(nstencil_type);
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = 0; j <= atom->ntypes; j++) {
      if (maxstencil_type[i][j] > 0) memory->destroy(stencil_type[i][j]);
    }
    delete [] stencil_type[i];
  }
  delete [] stencil_type;
  memory->destroy(maxstencil_type);
}

void NStencilHalfBytype3dNewtoff::copy_bin_info_bytype(int itype) {

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

int NStencilHalfBytype3dNewtoff::copy_neigh_info_bytype(int itype) {

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

void NStencilHalfBytype3dNewtoff::create_setup()
{

  int itype, jtype;
  int maxtypes;
  int smax;

  //printf("NStencilHalfBytype3dNewtoff::create_steup()\n");
  maxtypes = atom->ntypes;

  if (maxstencil_type == NULL) {
    memory->create(maxstencil_type, maxtypes+1, maxtypes+1, "BAD A");
    memory->create(nstencil_type, maxtypes+1, maxtypes+1, "BAD B");
    stencil_type = new int**[maxtypes+1]();
    for (itype = 1; itype <= maxtypes; ++itype) {
      stencil_type[itype] = new int*[maxtypes+1]();
      for (jtype = 1; jtype <= maxtypes; ++jtype) {
	maxstencil_type[itype][jtype] = 0;
      }
    }
  }

  // like -> like => use standard newtoff stencil at bin 

  for (itype = 1; itype <= maxtypes; ++itype) {
    copy_bin_info_bytype(itype);
    smax = copy_neigh_info_bytype(itype);
    if (smax > maxstencil_type[itype][itype]) {
      maxstencil_type[itype][itype] = smax;
      memory->destroy(stencil_type[itype][itype]);
      memory->create(stencil_type[itype][itype], smax,
		     "NStencilHalfBytypeNewtoff::create_steup() stencil");
    }
    create_newtoff(itype, itype, cutneighmaxsq);
  }

  // smaller -> larger => use existing newtoff stencil in larger bin
  // larger -> smaller => use multi-like stencil for small-large in smaller bin
  // If types are same cutoff, use existing like-like stencil.

  for (itype = 1; itype <= maxtypes; ++itype) {
    for (jtype = 1; jtype <= maxtypes; ++jtype) {
      if (itype == jtype) continue;
      if (cuttypesq[itype] <= cuttypesq[jtype]) {
	// Potential destroy/create problem?
	nstencil_type[itype][jtype] = nstencil_type[jtype][jtype];
	stencil_type[itype][jtype]  = stencil_type[jtype][jtype];
      }
      else {
	copy_bin_info_bytype(jtype);
	// smax = copy_neigh_info_bytype(jtype);

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
	  memory->create(stencil_type[itype][jtype], smax, "Bad C");
	}
	create_newtoff(itype, jtype, cuttypesq[jtype]);
      }
    }
  }

  //for (itype = 1; itype <= maxtypes; itype++) {
  //  for (jtype = 1; jtype <= maxtypes; jtype++) {
  //    printf("i j n %d %d %d\n", itype, jtype, nstencil_type[itype][jtype]);
  //    printf("i j n %d %d %d\n", itype, jtype, maxstencil_type[itype][jtype]);
  //  }
  // }

}

void NStencilHalfBytype3dNewtoff::create_newtoff(int itype, int jtype, double cutsq) {

  int i, j, k, ns;

  ns = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
	if (bin_distance(i,j,k) < cutsq)
	  stencil_type[itype][jtype][ns++] = k*mbiny*mbinx + j*mbinx + i;

  nstencil_type[itype][jtype] = ns;
}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilHalfBytype3dNewtoff::create()
{
  //int i,j,k;

  //nstencil = 0;

  //for (k = -sz; k <= sz; k++)
  //  for (j = -sy; j <= sy; j++)
  //    for (i = -sx; i <= sx; i++)
  //      if (bin_distance(i,j,k) < cutneighmaxsq)
  //        stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
}

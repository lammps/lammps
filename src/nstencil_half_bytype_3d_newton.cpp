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

#include "nstencil_half_bytype_3d_newton.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "nbin.h"
#include "memory.h"
#include "atom.h"

#include <math.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfBytype3dNewton::NStencilHalfBytype3dNewton(LAMMPS *lmp) :
  NStencil(lmp)
{
  maxstencil_type = NULL;
}

NStencilHalfBytype3dNewton::~NStencilHalfBytype3dNewton() {

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

// KS To superclass
void NStencilHalfBytype3dNewton::copy_bin_info_bytype(int itype) {

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

// KS To superclass?
int NStencilHalfBytype3dNewton::copy_neigh_info_bytype(int itype) {

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

void NStencilHalfBytype3dNewton::create_setup()
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

void NStencilHalfBytype3dNewton::create_newton(int itype, int jtype, double cutsq) {

  int i, j, k, ns;

  ns = 0;

  for (k = 0; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
	if (k > 0 || j > 0 || (j == 0 && i > 0)) { 
	  if (bin_distance(i,j,k) < cutsq)
	    stencil_type[itype][jtype][ns++] = k*mbiny*mbinx + j*mbinx + i;
	  }
  nstencil_type[itype][jtype] = ns;
}

void NStencilHalfBytype3dNewton::create_newtoff(int itype, int jtype, double cutsq) {

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

void NStencilHalfBytype3dNewton::create()
{
  // KS. Move "creation" here.
}

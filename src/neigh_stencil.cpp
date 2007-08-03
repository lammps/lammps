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
#include "memory.h"

using namespace LAMMPS_NS;

enum{NSQ,BIN,MULTI};     // also in neighbor.cpp

/* ----------------------------------------------------------------------
   routines to create a stencil = list of bin offsets
   stencil = bins whose closest corner to central bin is within cutoff
   sx,sy,sz = bin bounds = furthest the stencil could possibly extend
   3d creates xyz stencil, 2d creates xy stencil
   for half neigh list with partial Newton:
     stencil is all surrounding bins
     stencil includes self
   for half neigh list with full Newton:
     stencil is bins to the "upper right" of central bin
     stencil does not include self
   for half neigh list with triclinic:
     stencil is all bins in z-plane of self and above, but not below
     in 2d is all bins in y-plane of self and above, but not below
     stencil includes self
   for full neigh list:
     stencil is all surrounding bins including self
     regardless of newton on/off or triclinic
   for multi:
     create one stencil for each atom type
     stencil is same bin bounds as newton on/off, triclinic, half/full
     cutoff is not cutneighmaxsq, but max cutoff for that atom type
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_allocate(int sx, int sy, int sz)
{
  int i;

  // for 2d, sz = 0

  int nmax = (2*sz+1) * (2*sy+1) * (2*sx+1);

  if (half) {
    if (style == BIN) {
      if (nmax > maxstencil) {
	maxstencil = nmax;
	memory->sfree(stencil);
	stencil = (int *) memory->smalloc(nmax*sizeof(int),"neigh:stencil");
      }

    } else {
      int n = atom->ntypes;
      if (nstencil_multi == NULL) {
	maxstencil_multi = 0;
	nstencil_multi = new int[n+1];
	stencil_multi = new int*[n+1];
	distsq_multi = new double*[n+1];
	for (i = 1; i <= n; i++) {
	  nstencil_multi[i] = 0;
	  stencil_multi[i] = NULL;
	  distsq_multi[i] = NULL;
	}
      }
      if (nmax > maxstencil_multi) {
	maxstencil_multi = nmax;
	for (int i = 1; i <= n; i++) {
	  memory->sfree(stencil_multi[i]);
	  memory->sfree(distsq_multi[i]);
	  stencil_multi[i] = (int *)
	    memory->smalloc(nmax*sizeof(int),"neigh:stencil_multi");
	  distsq_multi[i] = (double *) memory->smalloc(nmax*sizeof(double),
						       "neigh:distsq_multi");
	}
      }
    }
  }

  if (full) {
    if (style == BIN) {
      if (nmax > maxstencil_full) {
	maxstencil_full = nmax;
	memory->sfree(stencil_full);
	stencil_full = (int *) memory->smalloc(nmax*sizeof(int),
					       "neigh:stencil_full");
      }

    } else {
      int n = atom->ntypes;
      if (nstencil_full_multi == NULL) {
	maxstencil_full_multi = 0;
	nstencil_full_multi = new int[n+1];
	stencil_full_multi = new int*[n+1];
	distsq_full_multi = new double*[n+1];
	for (i = 1; i <= n; i++) {
	  nstencil_full_multi[i] = 0;
	  stencil_full_multi[i] = NULL;
	  distsq_full_multi[i] = NULL;
	}
      }
      if (nmax > maxstencil_full_multi) {
	maxstencil_full_multi = nmax;
	for (int i = 1; i <= n; i++) {
	  memory->sfree(stencil_full_multi[i]);
	  memory->sfree(distsq_full_multi[i]);
	  stencil_full_multi[i] = (int *)
	    memory->smalloc(nmax*sizeof(int),"neigh:stencil_full_multi");
	  distsq_full_multi[i] = 
	    (double *) memory->smalloc(nmax*sizeof(double),
				       "neigh:distsq_full_multi");
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_none(int sx, int sy, int sz) {}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_3d_no_newton(int sx, int sy, int sz)
{
  int i,j,k;
  nstencil = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
	if (bin_distance(i,j,k) < cutneighmaxsq)
	  stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_3d_no_newton_multi(int sx, int sy, int sz)
{
  int i,j,k,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (k = -sz; k <= sz; k++)
      for (j = -sy; j <= sy; j++)
	for (i = -sx; i <= sx; i++) {
	  rsq = bin_distance(i,j,k);
	  if (rsq < typesq) {
	    distsq[n] = rsq;
	    s[n++] = k*mbiny*mbinx + j*mbinx + i;
	  }
	}
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_3d_newton(int sx, int sy, int sz)
{
  int i,j,k;
  nstencil = 0;

  for (k = 0; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
	if (k > 0 || j > 0 || (j == 0 && i > 0))
	  if (bin_distance(i,j,k) < cutneighmaxsq)
	    stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_3d_newton_multi(int sx, int sy, int sz)
{
  int i,j,k,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (k = 0; k <= sz; k++)
      for (j = -sy; j <= sy; j++)
	for (i = -sx; i <= sx; i++)
	  if (k > 0 || j > 0 || (j == 0 && i > 0)) {
	    rsq = bin_distance(i,j,k);
	    if (rsq < typesq) {
	      distsq[n] = rsq;
	      s[n++] = k*mbiny*mbinx + j*mbinx + i;
	    }
	  }
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_3d_newton_tri(int sx, int sy, int sz)
{
  int i,j,k;
  nstencil = 0;

  for (k = 0; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
	if (bin_distance(i,j,k) < cutneighmaxsq)
	  stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_3d_newton_multi_tri(int sx, int sy, int sz)
{
  int i,j,k,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (k = 0; k <= sz; k++)
      for (j = -sy; j <= sy; j++)
	for (i = -sx; i <= sx; i++) {
	  rsq = bin_distance(i,j,k);
	  if (rsq < typesq) {
	    distsq[n] = rsq;
	    s[n++] = k*mbiny*mbinx + j*mbinx + i;
	  }
	}
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_2d_no_newton(int sx, int sy, int sz)
{
  int i,j;
  nstencil = 0;

  for (j = -sy; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq)
	stencil[nstencil++] = j*mbinx + i;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_2d_no_newton_multi(int sx, int sy, int sz)
{
  int i,j,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++) {
	rsq = bin_distance(i,j,0);
	if (rsq < typesq) {
	  distsq[n] = rsq;
	  s[n++] = j*mbinx + i;
	}
      }
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_2d_newton(int sx, int sy, int sz)
{
  int i,j;
  nstencil = 0;

  for (j = 0; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (j > 0 || (j == 0 && i > 0))
	if (bin_distance(i,j,0) < cutneighmaxsq)
	  stencil[nstencil++] = j*mbinx + i;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_2d_newton_multi(int sx, int sy, int sz)
{
  int i,j,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (j = 0; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
	if (j > 0 || (j == 0 && i > 0)) {
	  rsq = bin_distance(i,j,0);
	  if (rsq < typesq) {
	    distsq[n] = rsq;
	    s[n++] = j*mbinx + i;
	  }
	}
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_2d_newton_tri(int sx, int sy, int sz)
{
  int i,j;
  nstencil = 0;

  for (j = 0; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq)
	stencil[nstencil++] = j*mbinx + i;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_half_2d_newton_multi_tri(int sx, int sy, int sz)
{
  int i,j,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (j = 0; j <= sy; j++)
      for (i = -sx; i <= sx; i++) {
	rsq = bin_distance(i,j,0);
	if (rsq < typesq) {
	  distsq[n] = rsq;
	  s[n++] = j*mbinx + i;
	}
      }
    nstencil_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_full_3d(int sx, int sy, int sz)
{
  int i,j,k;
  nstencil_full = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
	if (bin_distance(i,j,k) < cutneighmaxsq)
	  stencil_full[nstencil_full++] = k*mbiny*mbinx + j*mbinx + i;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_full_3d_multi(int sx, int sy, int sz)
{
  int i,j,k,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_full_multi[itype];
    distsq = distsq_full_multi[itype];
    n = 0;
    for (k = -sz; k <= sz; k++)
      for (j = -sy; j <= sy; j++)
	for (i = -sx; i <= sx; i++) {
	  rsq = bin_distance(i,j,k);
	  if (rsq < typesq) {
	    distsq[n] = rsq;
	    s[n++] = k*mbiny*mbinx + j*mbinx + i;
	  }
	}
    nstencil_full_multi[itype] = n;
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_full_2d(int sx, int sy, int sz)
{
  int i,j;
  nstencil_full = 0;

  for (j = -sy; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq)
	stencil_full[nstencil_full++] = j*mbinx + i;
}

/* ---------------------------------------------------------------------- */

void Neighbor::stencil_full_2d_multi(int sx, int sy, int sz)
{
  int i,j,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_full_multi[itype];
    distsq = distsq_full_multi[itype];
    n = 0;
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++) {
	rsq = bin_distance(i,j,0);
	if (rsq < typesq) {
	  distsq[n] = rsq;
	  s[n++] = j*mbinx + i;
	}
      }
    nstencil_full_multi[itype] = n;
  }
}

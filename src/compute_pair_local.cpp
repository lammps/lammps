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

#include "math.h"
#include "string.h"
#include "compute_pair_local.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputePairLocal::ComputePairLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal compute pair/local command");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  dflag = eflag = fflag = -1;
  nvalues = 0;

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;
    if (strcmp(arg[iarg],"dist") == 0) dflag = nvalues++;
    else if (strcmp(arg[iarg],"eng") == 0) eflag = nvalues++;
    else if (strcmp(arg[iarg],"force") == 0) fflag = nvalues++;
    else error->all("Invalid keyword in compute pair/local command");
  }

  nmax = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputePairLocal::~ComputePairLocal()
{
  memory->sfree(vector);
  memory->destroy_2d_double_array(array);
}

/* ---------------------------------------------------------------------- */

void ComputePairLocal::init()
{
  if (force->pair == NULL) 
    error->all("No pair style is defined for compute pair/local");
  if (force->pair->single_enable == 0)
    error->all("Pair style does not support compute pair/local");

  // need an occasional half neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputePairLocal::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePairLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute pair info

  ncount = compute_pairs(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_pairs(1);
}

/* ----------------------------------------------------------------------
   count pairs and compute pair info on this proc
   only count pair once if newton_pair is off
   both atom I,J must be in group
   if flag is set, compute requested info about pair
------------------------------------------------------------------------- */

int ComputePairLocal::compute_pairs(int flag)
{
  int i,j,m,n,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,eng,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double *dbuf,*ebuf,*fbuf;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)

  if (flag == 0) neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I or J are not in group

  if (flag) {
    if (nvalues == 1) {
      if (dflag >= 0) dbuf = vector;
      if (eflag >= 0) ebuf = vector;
      if (fflag >= 0) fbuf = vector;
    } else {
      if (dflag >= 0) dbuf = &array[0][dflag];
      if (eflag >= 0) ebuf = &array[0][eflag];
      if (fflag >= 0) fbuf = &array[0][fflag];
    }
  }

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  m = n = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      if (!(mask[j] & groupbit)) continue;
      if (newton_pair == 0 && j >= nlocal) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;
	
      if (flag) {
	if (dflag >= 0) dbuf[n] = sqrt(rsq);
	if (eflag >= 0 || fflag >= 0) {
	  eng = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
	  if (eflag >= 0) ebuf[n] = eng;
	  if (fflag >= 0) fbuf[n] = fpair;
	}
	n += nvalues;
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputePairLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->sfree(vector);
    vector = (double *) memory->smalloc(nmax*sizeof(double),
					"pair/local:vector");
    vector_local = vector;
  } else {
    memory->destroy_2d_double_array(array);
    array = memory->create_2d_double_array(nmax,nvalues,
					   "pair/local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePairLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}

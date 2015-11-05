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

/* ----------------------------------------------------------------------
   Contributing author:  Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <complex>
#include <string.h>
#include <stdlib.h>
#include "compute_hexorder_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeHexOrderAtom::ComputeHexOrderAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3 ) error->all(FLERR,"Illegal compute hexorder/atom command");

  nnn = 6;

  // process optional args

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"degree") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal lattice command");
      nnn = force->numeric(FLERR,arg[iarg+1]);
      if (nnn < 0)
        error->all(FLERR,"Illegal lattice command");
      iarg += 2;
    }
  }

  ncol = 2;
  peratom_flag = 1;
  size_peratom_cols = ncol;

  nmax = 0;
  q6array = NULL;
  maxneigh = 0;
  distsq = NULL;
  nearest = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeHexOrderAtom::~ComputeHexOrderAtom()
{
  memory->destroy(q6array);
  memory->destroy(distsq);
  memory->destroy(nearest);
}

/* ---------------------------------------------------------------------- */

void ComputeHexOrderAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute hexorder/atom requires a pair style be defined");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"hexorder/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute hexorder/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeHexOrderAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeHexOrderAtom::compute_peratom()
{
  int i,j,m,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow order parameter array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(q6array);
    nmax = atom->nmax;
    memory->create(q6array,nmax,ncol,"hexorder/atom:q6array");
    array_atom = q6array;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute order parameter for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *mask = atom->mask;
  double cutsq = force->pair->cutforce * force->pair->cutforce;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double* q6 = q6array[i];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      
      // insure distsq and nearest arrays are long enough

      if (jnum > maxneigh) {
        memory->destroy(distsq);
        memory->destroy(nearest);
        maxneigh = jnum;
        memory->create(distsq,maxneigh,"hexcoord/atom:distsq");
        memory->create(nearest,maxneigh,"hexcoord/atom:nearest");
      }

      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // nearest[] = atom indices of neighbors

      int ncount = 0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          distsq[ncount] = rsq;
          nearest[ncount++] = j;
        }
      }

      // if not nnn neighbors, order parameter = 0;

      if (ncount < nnn) {
	q6[0] = q6[1] = 0.0;
        continue;
      }

      // store nnn nearest neighs in 1st nnn locations of distsq and nearest

      select2(nnn,ncount,distsq,nearest);

      double usum = 0.0;
      double vsum = 0.0;
      
      for (jj = 0; jj < nnn; jj++) {
	j = nearest[jj];
	j &= NEIGHMASK;
	
	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	double u, v;
	calc_qn(delx, dely, u, v);
	usum += u;
	vsum += v;
      }
      q6[0] = usum/nnn;
      q6[1] = vsum/nnn;
    }
  }
}

// this might be faster than pow(std::complex) on some platforms

inline void ComputeHexOrderAtom::calc_q6(double delx, double dely, double &u, double &v) {
  double rinv = 1.0/sqrt(delx*delx+dely*dely);
  double x = delx*rinv;
  double y = dely*rinv;
  double a = x*x;
  double b1 = y*y;
  double b2 = b1*b1;
  double b3 = b2*b1;

  // (x + i y)^6 coeffs: 1, 6, -15, -20, 15, 6, -1

  u = ((  a - 15*b1)*a + 15*b2)*a - b3;
  v = ((6*a - 20*b1)*a +  6*b2)*x*y;
}

inline void ComputeHexOrderAtom::calc_qn(double delx, double dely, double &u, double &v) {
  double rinv = 1.0/sqrt(delx*delx+dely*dely);
  double x = delx*rinv;
  double y = dely*rinv;
  std::complex<double> z = x + y*1i;
  std::complex<double> zn = pow(z,nnn);
  u = real(zn);
  v = imag(zn);
}

/* ----------------------------------------------------------------------
   select2 routine from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   sort auxiliary array at same time
------------------------------------------------------------------------- */

#define SWAP(a,b)   tmp = a; a = b; b = tmp;
#define ISWAP(a,b) itmp = a; a = b; b = itmp;

/* ---------------------------------------------------------------------- */

void ComputeHexOrderAtom::select2(int k, int n, double *arr, int *iarr)
{
  int i,ir,j,l,mid,ia,itmp;
  double a,tmp;

  arr--;
  iarr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr[l],arr[ir])
        ISWAP(iarr[l],iarr[ir])
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
      ISWAP(iarr[mid],iarr[l+1])
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir])
        ISWAP(iarr[l],iarr[ir])
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir])
        ISWAP(iarr[l+1],iarr[ir])
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1])
        ISWAP(iarr[l],iarr[l+1])
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i],arr[j])
        ISWAP(iarr[i],iarr[j])
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      iarr[l+1] = iarr[j];
      iarr[j] = ia;
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeHexOrderAtom::memory_usage()
{
  double bytes = ncol*nmax * sizeof(double);
  return bytes;
}

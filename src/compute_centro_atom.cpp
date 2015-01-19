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
   Contributing author: Michel Perez (U Lyon) for non-fcc lattices
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "compute_centro_atom.h"
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

ComputeCentroAtom::ComputeCentroAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute centro/atom command");

  if (strcmp(arg[3],"fcc") == 0) nnn = 12;
  else if (strcmp(arg[3],"bcc") == 0) nnn = 8;
  else nnn = force->inumeric(FLERR,arg[3]);

  if (nnn <= 0 || nnn % 2)
    error->all(FLERR,"Illegal neighbor value for compute centro/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  centro = NULL;
  maxneigh = 0;
  distsq = NULL;
  nearest = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeCentroAtom::~ComputeCentroAtom()
{
  memory->destroy(centro);
  memory->destroy(distsq);
  memory->destroy(nearest);
}

/* ---------------------------------------------------------------------- */

void ComputeCentroAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute centro/atom requires a pair style be defined");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"centro/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute centro/atom");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeCentroAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCentroAtom::compute_peratom()
{
  int i,j,k,ii,jj,kk,n,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,value;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow centro array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(centro);
    nmax = atom->nmax;
    memory->create(centro,nmax,"centro/atom:centro");
    vector_atom = centro;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // npairs = number of unique pairs

  int nhalf = nnn/2;
  int npairs = nnn * (nnn-1) / 2;
  double *pairs = new double[npairs];

  // compute centro-symmetry parameter for each atom in group
  // use full neighbor list

  double **x = atom->x;
  int *mask = atom->mask;
  double cutsq = force->pair->cutforce * force->pair->cutforce;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
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
        memory->create(distsq,maxneigh,"centro/atom:distsq");
        memory->create(nearest,maxneigh,"centro/atom:nearest");
      }

      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // nearest[] = atom indices of neighbors

      n = 0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          distsq[n] = rsq;
          nearest[n++] = j;
        }
      }

      // if not nnn neighbors, centro = 0.0

      if (n < nnn) {
        centro[i] = 0.0;
        continue;
      }

      // store nnn nearest neighs in 1st nnn locations of distsq and nearest

      select2(nnn,n,distsq,nearest);

      // R = Ri + Rj for each of npairs i,j pairs among nnn neighbors
      // pairs = squared length of each R

      n = 0;
      for (j = 0; j < nnn; j++) {
        jj = nearest[j];
        for (k = j+1; k < nnn; k++) {
          kk = nearest[k];
          delx = x[jj][0] + x[kk][0] - 2.0*xtmp;
          dely = x[jj][1] + x[kk][1] - 2.0*ytmp;
          delz = x[jj][2] + x[kk][2] - 2.0*ztmp;
          pairs[n++] = delx*delx + dely*dely + delz*delz;
        }
      }

      // store nhalf smallest pair distances in 1st nhalf locations of pairs

      select(nhalf,npairs,pairs);

      // centrosymmetry = sum of nhalf smallest squared values

      value = 0.0;
      for (j = 0; j < nhalf; j++) value += pairs[j];
      centro[i] = value;
    } else centro[i] = 0.0;
  }

  delete [] pairs;
}

/* ----------------------------------------------------------------------
   2 select routines from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   2nd routine sorts auxiliary array at same time
------------------------------------------------------------------------- */

#define SWAP(a,b)   tmp = a; a = b; b = tmp;
#define ISWAP(a,b) itmp = a; a = b; b = itmp;

void ComputeCentroAtom::select(int k, int n, double *arr)
{
  int i,ir,j,l,mid;
  double a,tmp;

  arr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr[l],arr[ir])
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir])
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir])
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1])
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i],arr[j])
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeCentroAtom::select2(int k, int n, double *arr, int *iarr)
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

double ComputeCentroAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}

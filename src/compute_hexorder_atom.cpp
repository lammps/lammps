// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#include "compute_hexorder_atom.h"
#include <cmath>
#include <cstring>
#include <complex>
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
#include "math_const.h"

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeHexOrderAtom::ComputeHexOrderAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  distsq(nullptr), nearest(nullptr), qnarray(nullptr)
{
  if (narg < 3 ) error->all(FLERR,"Illegal compute hexorder/atom command");

  ndegree = 6;
  nnn = 6;
  cutsq = 0.0;

  // process optional args

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"degree") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute hexorder/atom command");
      ndegree = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (ndegree < 0)
        error->all(FLERR,"Illegal compute hexorder/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nnn") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute hexorder/atom command");
      if (strcmp(arg[iarg+1],"NULL") == 0)
        nnn = 0;
      else {
        nnn = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (nnn < 0)
          error->all(FLERR,"Illegal compute hexorder/atom command");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute hexorder/atom command");
      double cutoff = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (cutoff <= 0.0)
        error->all(FLERR,"Illegal compute hexorder/atom command");
      cutsq = cutoff*cutoff;
      iarg += 2;
    } else error->all(FLERR,"Illegal compute hexorder/atom command");
  }

  ncol = 2;
  peratom_flag = 1;
  size_peratom_cols = ncol;

  nmax = 0;
  maxneigh = 0;
}

/* ---------------------------------------------------------------------- */

ComputeHexOrderAtom::~ComputeHexOrderAtom()
{
  memory->destroy(qnarray);
  memory->destroy(distsq);
  memory->destroy(nearest);
}

/* ---------------------------------------------------------------------- */

void ComputeHexOrderAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute hexorder/atom requires a pair style be defined");
  if (cutsq == 0.0) cutsq = force->pair->cutforce * force->pair->cutforce;
  else if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,
               "Compute hexorder/atom cutoff is longer than pairwise cutoff");

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

void ComputeHexOrderAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeHexOrderAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow order parameter array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(qnarray);
    nmax = atom->nmax;
    memory->create(qnarray,nmax,ncol,"hexorder/atom:qnarray");
    array_atom = qnarray;
  }

  // invoke full neighbor list (will copy or build if necessary)
  // on the first step of a run, set preflag to one in neighbor->build_one(...)

  if (update->firststep == update->ntimestep) neighbor->build_one(list,1);
  else neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute order parameter for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *mask = atom->mask;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double* qn = qnarray[i];
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
        memory->create(distsq,maxneigh,"hexorder/atom:distsq");
        memory->create(nearest,maxneigh,"hexorder/atom:nearest");
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
        qn[0] = qn[1] = 0.0;
        continue;
      }

      // if nnn > 0, use only nearest nnn neighbors

      if (nnn > 0) {
        select2(nnn,ncount,distsq,nearest);
        ncount = nnn;
      }

      double usum = 0.0;
      double vsum = 0.0;

      for (jj = 0; jj < ncount; jj++) {
        j = nearest[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        double u, v;
        calc_qn_complex(delx, dely, u, v);
        usum += u;
        vsum += v;
      }
      qn[0] = usum/nnn;
      qn[1] = vsum/nnn;
    } else qn[0] = qn[1] = 0.0;
  }
}

// calculate order parameter using std::complex::pow function

inline void ComputeHexOrderAtom::calc_qn_complex(double delx, double dely, double &u, double &v) {
  double rinv = 1.0/sqrt(delx*delx+dely*dely);
  double x = delx*rinv;
  double y = dely*rinv;
  std::complex<double> z(x, y);
  std::complex<double> zn = pow(z, ndegree);
  u = real(zn);
  v = imag(zn);
}

// calculate order parameter using trig functions
// this is usually slower, but can be used if <complex> not available

inline void ComputeHexOrderAtom::calc_qn_trig(double delx, double dely, double &u, double &v) {
  double ntheta;
  if (fabs(delx) <= MY_EPSILON) {
    if (dely > 0.0) ntheta = ndegree * MY_PI / 2.0;
    else ntheta = ndegree * 3.0 * MY_PI / 2.0;
  } else ntheta = ndegree * atan(dely / delx);
  u = cos(ntheta);
  v = sin(ntheta);
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
  double bytes = (double)ncol*nmax * sizeof(double);
  bytes += (double)maxneigh * sizeof(double);
  bytes += (double)maxneigh * sizeof(int);

  return bytes;
}

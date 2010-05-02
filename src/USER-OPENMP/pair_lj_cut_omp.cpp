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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_cut_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairLJCutOMP::PairLJCutOMP(LAMMPS *lmp) : PairOMP(lmp)
{
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairLJCutOMP::~PairLJCutOMP()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_double_array(cut);
    memory->destroy_2d_double_array(epsilon);
    memory->destroy_2d_double_array(sigma);
    memory->destroy_2d_double_array(lj1);
    memory->destroy_2d_double_array(lj2);
    memory->destroy_2d_double_array(lj3);
    memory->destroy_2d_double_array(lj4);
    memory->destroy_2d_double_array(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) return eval<1,1,1>();
      else return eval<1,1,0>();
    } else {
      if (force->newton_pair) return eval<1,0,1>();
      else return eval<1,0,0>();
    }
  } else {
    if (force->newton_pair) return eval<0,0,1>();
    else return eval<0,0,0>();
  }
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR> 
void PairLJCutOMP::eval() 
{

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    const int nthreads = comm->nthreads;

    int i,j,ii,jj,inum,jnum,itype,jtype;
    double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
    double rsq,r2inv,r6inv,forcelj,factor_lj;
    int *ilist,*jlist,*numneigh,**firstneigh;

    evdwl = 0.0;

    const int nlocal = atom->nlocal;
    const int nall = nlocal + atom->nghost;

    double **x = atom->x;
    // each thread operates on its own copy of the forces array
    double **f = atom->f + nall*tid;
    int *type = atom->type;
    double *special_lj = force->special_lj;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  
    // loop over neighbors of my atoms
    // each thread works on a fixed chunk of atoms.
    // XXX: no load balancing. 
    // need to check how to estimate equal amounts of work. 
    // add per thread profiling to see how large a difference.
    const int iidelta = (nthreads > 1) ? 1 + inum/nthreads : inum;
    const int iifrom = tid*iidelta;
    const int iito   = ((iifrom + iidelta) > inum) ? inum : (iifrom + iidelta);
    for (ii = iifrom; ii < iito; ++ii) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];

	if (j < nall) factor_lj = 1.0;
	else {
	  factor_lj = special_lj[j/nall];
	  j %= nall;
	}

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	jtype = type[j];

	if (rsq < cutsq[itype][jtype]) {
	  r2inv = 1.0/rsq;
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  fpair = factor_lj*forcelj*r2inv;

	  f[i][0] += delx*fpair;
	  f[i][1] += dely*fpair;
	  f[i][2] += delz*fpair;
	  if (NEWTON_PAIR || j < nlocal) {
	    f[j][0] -= delx*fpair;
	    f[j][1] -= dely*fpair;
	    f[j][2] -= delz*fpair;
	  }

	  if (EFLAG) {
	    evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) 
	      - offset[itype][jtype];
	    evdwl *= factor_lj;
	  }

	  if (EVFLAG) ev_tally_thr(i,j,nlocal,NEWTON_PAIR,
				   evdwl,0.0,fpair,delx,dely,delz,tid);
	}
      }
    }

    // reduce per thread forces into global force array.
    // post a barrier to wait until all threads are done.
    // the reduction can be threaded as well.
#if defined(_OPENMP)
#pragma omp barrier
    double **fall = atom->f;
    const int idelta = (nthreads > 1) ? 1 + nall/nthreads : nall;
    const int ifrom = tid*idelta;
    const int ito   = ((ifrom + idelta) > nall) ? nall : (ifrom + idelta);
    for (int n = 1; n < nthreads; ++n) {
      const int toffs = n*nall;
      f = fall + toffs;
      for (int m = ifrom; m < ito; ++m) {
	fall[m][0] += f[m][0];
	f[m][0] = 0.0;
	fall[m][1] += f[m][1];
	f[m][1] = 0.0;
	fall[m][2] += f[m][2];
	f[m][2] = 0.0;
      }
    }
#endif
  }
  ev_reduce_thr();
  if (vflag_fdotr) virial_compute();
}

/* ---------------------------------------------------------------------- */

void PairLJCutOMP::compute_inner()
{
  if (force->newton_pair) return eval_inner<1>();
  else return eval_inner<0>();
}

template <int NEWTON_PAIR>
void PairLJCutOMP::eval_inner() 
{
#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    const int nthreads = comm->nthreads;

    int i,j,ii,jj,inum,jnum,itype,jtype;
    double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
    double rsq,r2inv,r6inv,forcelj,factor_lj,rsw;
    int *ilist,*jlist,*numneigh,**firstneigh;

    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;

    double **x = atom->x;
    double **f = atom->f + nall*tid;
    int *type = atom->type;
    double *special_lj = force->special_lj;
  
    inum = listinner->inum;
    ilist = listinner->ilist;
    numneigh = listinner->numneigh;
    firstneigh = listinner->firstneigh;

    double cut_out_on = cut_respa[0];
    double cut_out_off = cut_respa[1];
  
    double cut_out_diff = cut_out_off - cut_out_on;
    double cut_out_on_sq = cut_out_on*cut_out_on;
    double cut_out_off_sq = cut_out_off*cut_out_off;
  
    // loop over neighbors of my atoms
    // each thread works on a fixed chunk of atoms.
    // XXX: no load balancing! see compute() method
    const int iidelta = (nthreads > 1) ? 1 + inum/nthreads : inum;
    const int iifrom = tid*iidelta;
    const int iito   = ((iifrom + iidelta) > inum) ? inum : (iifrom + iidelta);
    for (ii = iifrom; ii < iito; ++ii) {

      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];
      
	if (j < nall) factor_lj = 1.0;
	else {
	  factor_lj = special_lj[j/nall];
	  j %= nall;
	}

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq < cut_out_off_sq) {
	  r2inv = 1.0/rsq;
	  r6inv = r2inv*r2inv*r2inv;
	  jtype = type[j];
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  fpair = factor_lj*forcelj*r2inv;
	  if (rsq > cut_out_on_sq) {
	    rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff; 
	    fpair *= 1.0 - rsw*rsw*(3.0 - 2.0*rsw);
	  }

	  f[i][0] += delx*fpair;
	  f[i][1] += dely*fpair;
	  f[i][2] += delz*fpair;
	  if (NEWTON_PAIR || j < nlocal) {
	    f[j][0] -= delx*fpair;
	    f[j][1] -= dely*fpair;
	    f[j][2] -= delz*fpair;
	  }
	}
      }
    }

    // reduce per thread forces into global force array.
    // post a barrier to wait until all threads are done.
    // the reduction can be threaded as well.
#if defined(_OPENMP)
#pragma omp barrier
    double **fall = atom->f;
    const int idelta = (nthreads > 1) ? 1 + nall/nthreads : nall;
    const int ifrom = tid*idelta;
    const int ito   = ((ifrom + idelta) > nall) ? nall : (ifrom + idelta);
    for (int n = 1; n < nthreads; ++n) {
      const int toffs = n*nall;
      f = fall + toffs;
      for (int m = ifrom; m < ito; ++m) {
	fall[m][0] += f[m][0];
	f[m][0] = 0.0;
	fall[m][1] += f[m][1];
	f[m][1] = 0.0;
	fall[m][2] += f[m][2];
	f[m][2] = 0.0;
      }
    }
#endif
  }
}

/* ---------------------------------------------------------------------- */


void PairLJCutOMP::compute_middle()
{
  if (force->newton_pair) return eval_middle<1>();
  else return eval_middle<0>();
}

template <int NEWTON_PAIR>
void PairLJCutOMP::eval_middle() 
{
#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    const int nthreads = comm->nthreads;

    int i,j,ii,jj,inum,jnum,itype,jtype;
    double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
    double rsq,r2inv,r6inv,forcelj,factor_lj,rsw;
    int *ilist,*jlist,*numneigh,**firstneigh;

    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;
    double **x = atom->x;
    double **f = atom->f + nall*tid;
    int *type = atom->type;
    double *special_lj = force->special_lj;

    inum = listmiddle->inum;
    ilist = listmiddle->ilist;
    numneigh = listmiddle->numneigh;
    firstneigh = listmiddle->firstneigh;

    double cut_in_off = cut_respa[0];
    double cut_in_on = cut_respa[1];
    double cut_out_on = cut_respa[2];
    double cut_out_off = cut_respa[3];

    double cut_in_diff = cut_in_on - cut_in_off;
    double cut_out_diff = cut_out_off - cut_out_on;
    double cut_in_off_sq = cut_in_off*cut_in_off;
    double cut_in_on_sq = cut_in_on*cut_in_on;
    double cut_out_on_sq = cut_out_on*cut_out_on;
    double cut_out_off_sq = cut_out_off*cut_out_off;

    // loop over neighbors of my atoms
    // each thread works on a fixed chunk of atoms.
    // XXX: no load balancing! see compute() method
    const int iidelta = (nthreads > 1) ? 1 + inum/nthreads : inum;
    const int iifrom = tid*iidelta;
    const int iito   = ((iifrom + iidelta) > inum) ? inum : (iifrom + iidelta);

    for (ii = iifrom; ii < iito; ++ii) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];

	if (j < nall) factor_lj = 1.0;
	else {
	  factor_lj = special_lj[j/nall];
	  j %= nall;
	}

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq < cut_out_off_sq && rsq > cut_in_off_sq) {
	  r2inv = 1.0/rsq;
	  r6inv = r2inv*r2inv*r2inv;
	  jtype = type[j];
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  fpair = factor_lj*forcelj*r2inv;
	  if (rsq < cut_in_on_sq) {
	    rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff; 
	    fpair *= rsw*rsw*(3.0 - 2.0*rsw);
	  }
	  if (rsq > cut_out_on_sq) {
	    rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff; 
	    fpair *= 1.0 + rsw*rsw*(2.0*rsw - 3.0);
	  }

	  f[i][0] += delx*fpair;
	  f[i][1] += dely*fpair;
	  f[i][2] += delz*fpair;
	  if (NEWTON_PAIR || j < nlocal) {
	    f[j][0] -= delx*fpair;
	    f[j][1] -= dely*fpair;
	    f[j][2] -= delz*fpair;
	  }
	}
      }
    }
    
    // reduce per thread forces into global force array.
    // post a barrier to wait until all threads are done.
    // the reduction can be threaded as well.
#if defined(_OPENMP)
#pragma omp barrier
    double **fall = atom->f;
    const int idelta = (nthreads > 1) ? 1 + nall/nthreads : nall;
    const int ifrom = tid*idelta;
    const int ito   = ((ifrom + idelta) > nall) ? nall : (ifrom + idelta);
    for (int n = 1; n < nthreads; ++n) {
      const int toffs = n*nall;
      f = fall + toffs;
      for (int m = ifrom; m < ito; ++m) {
	fall[m][0] += f[m][0];
	f[m][0] = 0.0;
	fall[m][1] += f[m][1];
	f[m][1] = 0.0;
	fall[m][2] += f[m][2];
	f[m][2] = 0.0;
      }
    }
#endif
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutOMP::compute_outer(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  if (evflag) {
    if (eflag) {
      if (vflag) {
	if (force->newton_pair) return eval_outer<1,1,1,1>();
	else return eval_outer<1,1,1,0>();
      } else {
	if (force->newton_pair) return eval_outer<1,1,0,1>();
	else return eval_outer<1,1,0,0>();
      }
    } else {
      if (vflag) {
	if (force->newton_pair) return eval_outer<1,0,1,1>();
	else return eval_outer<1,0,1,0>();
      } else {
	if (force->newton_pair) return eval_outer<1,0,0,1>();
	else return eval_outer<1,0,0,0>();
      }
    }
  } else {
    if (force->newton_pair) return eval_outer<0,0,0,1>();
    else return eval_outer<0,0,0,0>();
  }
}

template <int EVFLAG, int EFLAG, int VFLAG, int NEWTON_PAIR> 
void PairLJCutOMP::eval_outer() 
{
#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    const int nthreads = comm->nthreads;

    int i,j,ii,jj,inum,jnum,itype,jtype;
    double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
    double rsq,r2inv,r6inv,forcelj,factor_lj,rsw;
    int *ilist,*jlist,*numneigh,**firstneigh;

    evdwl = 0.0;

    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;
    double *special_lj = force->special_lj;

    inum = listouter->inum;
    ilist = listouter->ilist;
    numneigh = listouter->numneigh;
    firstneigh = listouter->firstneigh;

    double cut_in_off = cut_respa[2];
    double cut_in_on = cut_respa[3];

    double cut_in_diff = cut_in_on - cut_in_off;
    double cut_in_off_sq = cut_in_off*cut_in_off;
    double cut_in_on_sq = cut_in_on*cut_in_on;

    // loop over neighbors of my atoms
    // each thread works on a fixed chunk of atoms.
    // XXX: no load balancing! see compute() method
    const int iidelta = (nthreads > 1) ? 1 + inum/nthreads : inum;
    const int iifrom = tid*iidelta;
    const int iito   = ((iifrom + iidelta) > inum) ? inum : (iifrom + iidelta);

    for (ii = iifrom; ii < iito; ++ii) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];

	if (j < nall) factor_lj = 1.0;
	else {
	  factor_lj = special_lj[j/nall];
	  j %= nall;
	}

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	jtype = type[j];

	if (rsq < cutsq[itype][jtype]) {
	  if (rsq > cut_in_off_sq) {
	    r2inv = 1.0/rsq;
	    r6inv = r2inv*r2inv*r2inv;
	    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	    fpair = factor_lj*forcelj*r2inv;
	    if (rsq < cut_in_on_sq) {
	      rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff; 
	      fpair *= rsw*rsw*(3.0 - 2.0*rsw);
	    }

	    f[i][0] += delx*fpair;
	    f[i][1] += dely*fpair;
	    f[i][2] += delz*fpair;
	    if (NEWTON_PAIR || j < nlocal) {
	      f[j][0] -= delx*fpair;
	      f[j][1] -= dely*fpair;
	      f[j][2] -= delz*fpair;
	    }
	  }

	  if (EFLAG) {
	    r2inv = 1.0/rsq;
	    r6inv = r2inv*r2inv*r2inv;
	    evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
	      offset[itype][jtype];
	    evdwl *= factor_lj;
	  }

	  if (VFLAG) {
	    if (rsq <= cut_in_off_sq) {
	      r2inv = 1.0/rsq;
	      r6inv = r2inv*r2inv*r2inv;
	      forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	      fpair = factor_lj*forcelj*r2inv;
	    } else if (rsq < cut_in_on_sq)
	      fpair = factor_lj*forcelj*r2inv;
	  }

	  if (EVFLAG) ev_tally_thr(i,j,nlocal,NEWTON_PAIR,
				   evdwl,0.0,fpair,delx,dely,delz,0);
	}
      }
    }
    // reduce per thread forces into global force array.
    // post a barrier to wait until all threads are done.
    // the reduction can be threaded as well.
#if defined(_OPENMP)
#pragma omp barrier
    double **fall = atom->f;
    const int idelta = (nthreads > 1) ? 1 + nall/nthreads : nall;
    const int ifrom = tid*idelta;
    const int ito   = ((ifrom + idelta) > nall) ? nall : (ifrom + idelta);
    for (int n = 1; n < nthreads; ++n) {
      const int toffs = n*nall;
      f = fall + toffs;
      for (int m = ifrom; m < ito; ++m) {
	fall[m][0] += f[m][0];
	f[m][0] = 0.0;
	fall[m][1] += f[m][1];
	f[m][1] = 0.0;
	fall[m][2] += f[m][2];
	f[m][2] = 0.0;
      }
    }
#endif
  }
  ev_reduce_thr();
  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairLJCutOMP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  cut = memory->create_2d_double_array(n+1,n+1,"pair:cut");
  epsilon = memory->create_2d_double_array(n+1,n+1,"pair:epsilon");
  sigma = memory->create_2d_double_array(n+1,n+1,"pair:sigma");
  lj1 = memory->create_2d_double_array(n+1,n+1,"pair:lj1");
  lj2 = memory->create_2d_double_array(n+1,n+1,"pair:lj2");
  lj3 = memory->create_2d_double_array(n+1,n+1,"pair:lj3");
  lj4 = memory->create_2d_double_array(n+1,n+1,"pair:lj4");
  offset = memory->create_2d_double_array(n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairLJCutOMP::settings(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal pair_style command");

  cut_global = force->numeric(arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCutOMP::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(arg[2]);
  double sigma_one = force->numeric(arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = force->numeric(arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutOMP::init_style()
{
  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 1 && strcmp(update->integrate_style,"respa") == 0) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0) irequest = neighbor->request(this);
    else if (respa == 1) {
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    } else {
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 2;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respamiddle = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    }

  } else irequest = neighbor->request(this);

  // set rRESPA cutoffs

  if (strcmp(update->integrate_style,"respa") == 0 &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairLJCutOMP::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutOMP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
			       sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && cut[i][j] < cut_respa[3])
    error->all("Pair cutoff < Respa interior cutoff");

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);
        
    double PI = 4.0*atan(1.0);
    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*PI*all[0]*all[1]*epsilon[i][j] * 
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9); 
    ptail_ij = 16.0*PI*all[0]*all[1]*epsilon[i][j] * 
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9); 
  } 

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file 
------------------------------------------------------------------------- */

void PairLJCutOMP::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&epsilon[i][j],sizeof(double),1,fp);
	fwrite(&sigma[i][j],sizeof(double),1,fp);
	fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutOMP::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&epsilon[i][j],sizeof(double),1,fp);
	  fread(&sigma[i][j],sizeof(double),1,fp);
	  fread(&cut[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutOMP::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutOMP::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairLJCutOMP::single(int i, int j, int itype, int jtype, double rsq,
			    double factor_coul, double factor_lj,
			    double &fforce)
{
  double r2inv,r6inv,forcelj,philj;

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  fforce = factor_lj*forcelj*r2inv;

  philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
    offset[itype][jtype];
  return factor_lj*philj;
}

/* ---------------------------------------------------------------------- */

double PairLJCutOMP::memory_usage()
{
  const int n=atom->ntypes;
  
  double bytes = PairOMP::memory_usage();

  bytes += 9*((n+1)*(n+1) * sizeof(double) + (n+1)*sizeof(double *));
  bytes += 1*((n+1)*(n+1) * sizeof(int) + (n+1)*sizeof(int *));

  return bytes;
}

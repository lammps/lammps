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
#include "pair_lj_charmm_coul_charmm_implicit_omp.h"
#include "atom.h"
#include "force.h"
#include "neigh_list.h"
#include "comm.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulCharmmImplicitOMP::PairLJCharmmCoulCharmmImplicitOMP(LAMMPS *lmp) :
  PairLJCharmmCoulCharmmOMP(lmp)
{
  implicit = 1;
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmImplicitOMP::compute(int eflag, int vflag)
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
void PairLJCharmmCoulCharmmImplicitOMP::eval()
{

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int i,j,ii,jj,inum,jnum,itype,jtype,tid;
    double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
    double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
    double philj,switch1,switch2;
    int *ilist,*jlist,*numneigh,**firstneigh;

    evdwl = ecoul = 0.0;

    const int nlocal = atom->nlocal;
    const int nall = nlocal + atom->nghost;
    const int nthreads = comm->nthreads;

    double **x = atom->x;
    double **f = atom->f;
    double *q = atom->q;
    int *type = atom->type;
    double *special_coul = force->special_coul;
    double *special_lj = force->special_lj;
    double qqrd2e = force->qqrd2e;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    // loop over neighbors of my atoms

    int iifrom, iito;
    f = loop_setup_thr(f, iifrom, iito, tid, inum, nall, nthreads);
    for (ii = iifrom; ii < iito; ++ii) {
      i = ilist[ii];
      qtmp = q[i];
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

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq < cut_bothsq) {
	  r2inv = 1.0/rsq;

	  if (rsq < cut_coulsq) {
	    forcecoul = 2.0 * qqrd2e * qtmp*q[j]*r2inv;
	    if (rsq > cut_coul_innersq) {
	      switch1 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
		(cut_coulsq + 2.0*rsq - 3.0*cut_coul_innersq) / denom_coul;
	      switch2 = 12.0*rsq * (cut_coulsq-rsq) *
		(rsq-cut_coul_innersq) / denom_coul;
	      forcecoul *= switch1 + switch2;
	    }
	  } else forcecoul = 0.0;

	  if (rsq < cut_ljsq) {
	    r6inv = r2inv*r2inv*r2inv;
	    jtype = type[j];
	    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	    if (rsq > cut_lj_innersq) {
	      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		(cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	      switch2 = 12.0*rsq * (cut_ljsq-rsq) *
		(rsq-cut_lj_innersq) / denom_lj;
	      philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
	      forcelj = forcelj*switch1 + philj*switch2;
	    }
	  } else forcelj = 0.0;

	  fpair = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

	  f[i][0] += delx*fpair;
	  f[i][1] += dely*fpair;
	  f[i][2] += delz*fpair;
	  if (NEWTON_PAIR || j < nlocal) {
	    f[j][0] -= delx*fpair;
	    f[j][1] -= dely*fpair;
	    f[j][2] -= delz*fpair;
	  }

	  if (EFLAG) {
	    if (rsq < cut_coulsq) {
	      ecoul = qqrd2e * qtmp*q[j]*r2inv;
	      if (rsq > cut_coul_innersq) {
		switch1 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
		  (cut_coulsq + 2.0*rsq - 3.0*cut_coul_innersq) /
		  denom_coul;
		ecoul *= switch1;
	      }
	      ecoul *= factor_coul;
	    } else ecoul = 0.0;
	    if (rsq < cut_ljsq) {
	      evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
	      if (rsq > cut_lj_innersq) {
		switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		  (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
		evdwl *= switch1;
	      }
	      evdwl *= factor_lj;
	    } else evdwl = 0.0;
	  }

	  if (EVFLAG) ev_tally_thr(i,j,nlocal,NEWTON_PAIR,
				   evdwl,ecoul,fpair,delx,dely,delz,tid);
	}
      }
    }

    // reduce per thread forces and torques into global force/torque arrays.
    force_reduce_thr(atom->f, nall, nthreads, tid);
  }
  if (EVFLAG) ev_reduce_thr();
  if (vflag_fdotr) virial_compute();
}




/* ---------------------------------------------------------------------- */

double PairLJCharmmCoulCharmmImplicitOMP::single(int i, int j,
					      int itype, int jtype,
					      double rsq,
					      double factor_coul,
					      double factor_lj,
					      double &fforce)
{
  double r2inv,r6inv,switch1,switch2,forcecoul,forcelj,phicoul,philj;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq) {
    forcecoul = 2.0 * force->qqrd2e * atom->q[i]*atom->q[j]*r2inv;
    if (rsq > cut_coul_innersq) {
      switch1 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
	(cut_coulsq + 2.0*rsq - 3.0*cut_coul_innersq) / denom_coul;
      switch2 = 12.0*rsq * (cut_coulsq-rsq) * 
	(rsq-cut_coul_innersq) / denom_coul;
      forcecoul *= switch1 + switch2;
    }
  } else forcecoul = 0.0;
  if (rsq < cut_ljsq) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
    if (rsq > cut_lj_innersq) {
      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
	(cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
      switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
	(rsq-cut_lj_innersq) / denom_lj;
      philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
      forcelj = forcelj*switch1 + philj*switch2;
    }
  } else forcelj = 0.0;
  fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    phicoul = force->qqrd2e * atom->q[i]*atom->q[j]*r2inv;
    if (rsq > cut_coul_innersq) {
      switch1 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
	(cut_coulsq + 2.0*rsq - 3.0*cut_coul_innersq) / 
	denom_coul;
      phicoul *= switch1;
    }
    eng += factor_coul*phicoul;
  }
  if (rsq < cut_ljsq) {
    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
    if (rsq > cut_lj_innersq) {
      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
	(cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
      philj *= switch1;
    }
    eng += factor_lj*philj;
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairLJCharmmCoulCharmmImplicitOMP::extract(char *str, int &dim)
{
  void *rv = PairOMP::extract(str,dim);
  if (rv) return rv;

  dim = 2;
  if (strcmp(str,"lj14_1") == 0) return (void *) lj14_1;
  if (strcmp(str,"lj14_2") == 0) return (void *) lj14_2;
  if (strcmp(str,"lj14_3") == 0) return (void *) lj14_3;
  if (strcmp(str,"lj14_4") == 0) return (void *) lj14_4;

  dim = 0;
  if (strcmp(str,"implicit") == 0) return (void *) &implicit;
  return NULL;
}


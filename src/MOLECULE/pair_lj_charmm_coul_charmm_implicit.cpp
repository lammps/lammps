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
#include "pair_lj_charmm_coul_charmm_implicit.h"
#include "atom.h"
#include "force.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulCharmmImplicit::PairLJCharmmCoulCharmmImplicit(LAMMPS *lmp) :
  PairLJCharmmCoulCharmm(lmp) {}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmImplicit::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r2inv,r6inv,forcecoul,forcelj,fforce,factor_coul,factor_lj;
  double factor,phicoul,philj,switch1,switch2;
  int *ilist,*jlist,*numneigh,**firstneigh;

  eng_vdwl = eng_coul = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
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

	fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

	f[i][0] += delx*fforce;
	f[i][1] += dely*fforce;
	f[i][2] += delz*fforce;
	if (newton_pair || j < nlocal) {
	  f[j][0] -= delx*fforce;
	  f[j][1] -= dely*fforce;
	  f[j][2] -= delz*fforce;
	}

	if (eflag) {
	  if (newton_pair || j < nlocal) factor = 1.0;
	  else factor = 0.5;
	  if (rsq < cut_coulsq) {
	    phicoul = qqrd2e * qtmp*q[j]*r2inv;
	    if (rsq > cut_coul_innersq) {
	      switch1 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
		(cut_coulsq + 2.0*rsq - 3.0*cut_coul_innersq) / 
		denom_coul;
	      phicoul *= switch1;
	    }
	    eng_coul += factor*factor_coul*phicoul;
	  }
	  if (rsq < cut_ljsq) {
	    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
	    if (rsq > cut_lj_innersq) {
	      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		(cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	      philj *= switch1;
	    }
	    eng_vdwl += factor*factor_lj*philj;
	  }
	}

	if (vflag == 1) {
	  if (newton_pair == 0 && j >= nlocal) fforce *= 0.5;
	  virial[0] += delx*delx*fforce;
	  virial[1] += dely*dely*fforce;
	  virial[2] += delz*delz*fforce;
	  virial[3] += delx*dely*fforce;
	  virial[4] += delx*delz*fforce;
	  virial[5] += dely*delz*fforce;
	}
      }
    }
  }
  if (vflag == 2) virial_compute();
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmImplicit::single(int i, int j, int itype, int jtype,
					    double rsq, double factor_coul,
					    double factor_lj,
					    int eflag, One &one)
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
  one.fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;
  
  if (eflag) {
    if (rsq < cut_coulsq) {
      phicoul = force->qqrd2e * atom->q[i]*atom->q[j]*r2inv;
      if (rsq > cut_coul_innersq) {
	switch1 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
	  (cut_coulsq + 2.0*rsq - 3.0*cut_coul_innersq) / 
	  denom_coul;
	phicoul *= switch1;
      }
      one.eng_coul = factor_coul*phicoul;
    } else one.eng_coul = 0.0;
    if (rsq < cut_ljsq) {
      philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
      if (rsq > cut_lj_innersq) {
	switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
	  (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
	philj *= switch1;
      }
      one.eng_vdwl = factor_lj*philj;
    } else one.eng_vdwl = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmImplicit::extract_charmm(double ***p_lj14_1, 
						    double ***p_lj14_2,
						    double ***p_lj14_3,
						    double ***p_lj14_4,
						    int *p_implicit_flag)
{
  *p_lj14_1 = lj14_1;
  *p_lj14_2 = lj14_2;
  *p_lj14_3 = lj14_3;
  *p_lj14_4 = lj14_4;
  *p_implicit_flag = 1;
}

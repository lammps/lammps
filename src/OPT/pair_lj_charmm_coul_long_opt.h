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
   Contributing authors:  
     James Fischer, High Performance Technologies, Inc.
     David Richie, Stone Ridge Technology
     Vincent Natol, Stone Ridge Technology
------------------------------------------------------------------------- */

#ifndef PAIR_LJ_CHARMM_COUL_LONG_OPT_H
#define PAIR_LJ_CHARMM_COUL_LONG_OPT_H

#include "math.h"
#include "stdlib.h"
#include "pair_lj_charmm_coul_long.h"
#include "atom.h"
#include "force.h"
#include "neigh_list.h"

namespace LAMMPS_NS {

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define EWALD_A1  0.254829592
#define EWALD_A2 -0.284496736
#define EWALD_A3  1.421413741
#define EWALD_A4 -1.453152027
#define EWALD_A5  1.061405429

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

class PairLJCharmmCoulLongOpt : public PairLJCharmmCoulLong {
 public:
  PairLJCharmmCoulLongOpt(class LAMMPS *);
  void compute(int, int);

 private:
  template < int EVFLAG, int EFLAG, int NEWTON_PAIR > void eval();
};

template < int EVFLAG, int EFLAG, int NEWTON_PAIR >
void PairLJCharmmCoulLongOpt::eval()
{
  typedef struct { double x,y,z; } vec3_t;
  
  typedef struct {
    double cutsq,lj1,lj2,lj3,lj4,offset;
    double _pad[2];
  } fast_alpha_t;
  
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  double fraction,table;
  double r,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double philj,switch1,switch2;
 
  double rsq; 
  
  double evdwl = 0.0;
  double ecoul = 0.0;
  
  double** __restrict__ x = atom->x;
  double** __restrict__ f = atom->f;
  double* __restrict__ q = atom->q;
  int* __restrict__ type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double* __restrict__ special_coul = force->special_coul;
  double* __restrict__ special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  int* __restrict__ ilist = list->ilist;
  int** __restrict__ firstneigh = list->firstneigh;
  int* __restrict__ numneigh = list->numneigh;
  
  vec3_t* __restrict__ xx = (vec3_t*)x[0];
  vec3_t* __restrict__ ff = (vec3_t*)f[0];
  
  int ntypes = atom->ntypes;
  int ntypes2 = ntypes*ntypes;
  
  double tmp_coef1 = 1.0/denom_lj;
  double tmp_coef2 = cut_ljsq - 3.0*cut_lj_innersq;
  
  fast_alpha_t* __restrict__ fast_alpha = 
    (fast_alpha_t*)malloc(ntypes2*sizeof(fast_alpha_t));
  for (i = 0; i < ntypes; i++) for (j = 0; j < ntypes; j++) {
    fast_alpha_t& a = fast_alpha[i*ntypes+j];
    a.cutsq = cutsq[i+1][j+1];
    a.lj1 = lj1[i+1][j+1];
    a.lj2 = lj2[i+1][j+1];
    a.lj3 = lj3[i+1][j+1];
    a.lj4 = lj4[i+1][j+1];
  }
  fast_alpha_t* __restrict__ tabsix = fast_alpha;
  
  // loop over neighbors of my atoms
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double qtmp = q[i];
    double xtmp = xx[i].x;
    double ytmp = xx[i].y;
    double ztmp = xx[i].z;
    itype = type[i] - 1;
    int* __restrict__ jlist = firstneigh[i];
    jnum = numneigh[i];
    
    double tmpfx = 0.0;
    double tmpfy = 0.0;
    double tmpfz = 0.0;
    
    fast_alpha_t* __restrict__ tabsixi = (fast_alpha_t*) &tabsix[itype*ntypes];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      
      if (j < nall) {
	double delx = xtmp - xx[j].x;
	double dely = ytmp - xx[j].y;
	double delz = ztmp - xx[j].z;
	rsq = delx*delx + dely*dely + delz*delz;
	double tmp_coef3 = qtmp*q[j];
	
	if (rsq < cut_bothsq) {
	  r2inv = 1.0/rsq;
	  
	  forcecoul = 0.0;
	  if (rsq < cut_coulsq) {
	    if (!ncoultablebits || rsq <= tabinnersq) {
	      r = sqrt(rsq);
	      grij = g_ewald * r;
	      expm2 = exp(-grij*grij);
	      t = 1.0 / (1.0 + EWALD_P*grij);
	      erfc = t * 
		(EWALD_A1+t*(EWALD_A2+t*(EWALD_A3+t*(EWALD_A4+t*EWALD_A5)))) * 
		expm2;
	      prefactor = qqrd2e * tmp_coef3/r;
	      forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	    } else {
	      union_int_float_t rsq_lookup;
	      rsq_lookup.f = rsq;
	      itable = rsq_lookup.i & ncoulmask;
	      itable >>= ncoulshiftbits;
	      fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
	      table = ftable[itable] + fraction*dftable[itable];
	      forcecoul = tmp_coef3 * table;
	    }
	  }  
	  
	  forcelj = 0.0;
	  if (rsq < cut_ljsq) {
	    r6inv = r2inv*r2inv*r2inv;
	    jtype = type[j] - 1;
	    fast_alpha_t& a = tabsixi[jtype];
	    forcelj = r6inv * (a.lj1*r6inv - a.lj2);
	    if (rsq > cut_lj_innersq) {
	      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		(tmp_coef2 + 2.0*rsq) * tmp_coef1;
	      switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
		(rsq-cut_lj_innersq) * tmp_coef1;
	      philj = r6inv * (a.lj3*r6inv - a.lj4);
	      forcelj = forcelj*switch1 + philj*switch2;
	    }
	  } 
	  
	  double fpair = (forcecoul + forcelj) * r2inv;
	  
	  tmpfx += delx*fpair;
	  tmpfy += dely*fpair;
	  tmpfz += delz*fpair;
	  if (NEWTON_PAIR || j < nlocal) {
	    ff[j].x -= delx*fpair;
	    ff[j].y -= dely*fpair;
	    ff[j].z -= delz*fpair;
	  }
	  
	  if (EFLAG) {
	    if (rsq < cut_coulsq) {
	      if (!ncoultablebits || rsq <= tabinnersq)
		ecoul = prefactor*erfc;
	      else {
		table = etable[itable] + fraction*detable[itable];
		ecoul = tmp_coef3 * table;
	      }
	    } else ecoul = 0.0;

	    if (rsq < cut_ljsq) {
	      fast_alpha_t& a = tabsixi[jtype];
	      evdwl = r6inv*(a.lj3*r6inv-a.lj4);
	      if (rsq > cut_lj_innersq) {
		switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		  (tmp_coef2 + 2.0*rsq) * tmp_coef1;
		evdwl *= switch1;
	      }
	    } else evdwl = 0.0;
	  }

	  if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
			       evdwl,ecoul,fpair,delx,dely,delz);
	}

      } else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
	double delx = xtmp - xx[j].x;
	double dely = ytmp - xx[j].y;
	double delz = ztmp - xx[j].z;
	rsq = delx*delx + dely*dely + delz*delz;
	double tmp_coef3 = qtmp*q[j];
	
	if (rsq < cut_bothsq) {
	  r2inv = 1.0/rsq;
	  
	  forcecoul = 0.0;
	  if (rsq < cut_coulsq) {
	    if (!ncoultablebits || rsq <= tabinnersq) {
	      r = sqrt(rsq);
	      grij = g_ewald * r;
	      expm2 = exp(-grij*grij);
	      t = 1.0 / (1.0 + EWALD_P*grij);
	      erfc = t * 
		(EWALD_A1+t*(EWALD_A2+t*(EWALD_A3+t*(EWALD_A4+t*EWALD_A5)))) * 
		expm2;
	      prefactor = qqrd2e * tmp_coef3/r;
	      forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	      if (factor_coul < 1.0) {
		forcecoul -= (1.0-factor_coul)*prefactor;
	      }
	    } else {
	      union_int_float_t rsq_lookup;
	      rsq_lookup.f = rsq;
	      itable = rsq_lookup.i & ncoulmask;
	      itable >>= ncoulshiftbits;
	      fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
	      table = ftable[itable] + fraction*dftable[itable];
	      forcecoul = tmp_coef3 * table;
	      if (factor_coul < 1.0) {
		table = ctable[itable] + fraction*dctable[itable];
		prefactor = tmp_coef3 * table;
		forcecoul -= (1.0-factor_coul)*prefactor;
	      }
	    }
	  }  
	  
	  forcelj = 0.0;
	  if (rsq < cut_ljsq) {
	    r6inv = r2inv*r2inv*r2inv;
	    jtype = type[j] - 1;
	    fast_alpha_t& a = tabsixi[jtype];
	    forcelj = r6inv * (a.lj1*r6inv - a.lj2);
	    if (rsq > cut_lj_innersq) {
	      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		(tmp_coef2 + 2.0*rsq) * tmp_coef1;
	      switch2 = 12.0*rsq * (cut_ljsq-rsq) * 
		(rsq-cut_lj_innersq) * tmp_coef1;
	      fast_alpha_t& a = tabsixi[jtype];
	      philj = r6inv * (a.lj3*r6inv - a.lj4);
	      forcelj = forcelj*switch1 + philj*switch2;
	    }
	  }  
	  
	  double fpair = (forcecoul + factor_lj*forcelj) * r2inv;
	  
	  tmpfx += delx*fpair;
	  tmpfy += dely*fpair;
	  tmpfz += delz*fpair;
	  if (NEWTON_PAIR || j < nlocal) {
	    ff[j].x -= delx*fpair;
	    ff[j].y -= dely*fpair;
	    ff[j].z -= delz*fpair;
	  }
	  
	  if (EFLAG) {
	    if (rsq < cut_coulsq) {
	      if (!ncoultablebits || rsq <= tabinnersq)
		ecoul = prefactor*erfc;
	      else {
		table = etable[itable] + fraction*detable[itable];
		ecoul = tmp_coef3 * table;
	      }
	      if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
	    } else ecoul = 0.0;

	    if (rsq < cut_ljsq) {
	      fast_alpha_t& a = tabsixi[jtype];
	      evdwl = r6inv*(a.lj3*r6inv-a.lj4);
	      if (rsq > cut_lj_innersq) {
		switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
		  (tmp_coef2 + 2.0*rsq) * tmp_coef1;
		evdwl *= switch1;
	      }
	      evdwl *= factor_lj;
	    } else evdwl = 0.0;
	  }
	  
	  if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
			       evdwl,ecoul,fpair,delx,dely,delz);
	}
      }
    }

    ff[i].x += tmpfx;
    ff[i].y += tmpfy;
    ff[i].z += tmpfz;
  }

  free(fast_alpha); fast_alpha = 0;

  if (vflag_fdotr) virial_compute();
}

#undef EWALD_F
#undef EWALD_P
#undef EWALD_A1
#undef EWALD_A2
#undef EWALD_A3
#undef EWALD_A4
#undef EWALD_A5

#undef MIN
#undef MAX
 
}

#endif

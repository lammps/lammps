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
   Contributing authors: Randy Schunk (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "pair_yukawa_colloid_omp.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairYukawaColloidOMP::PairYukawaColloidOMP(LAMMPS *lmp) : PairYukawaOMP(lmp) {}

/* ---------------------------------------------------------------------- */

void PairYukawaColloidOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) return eval_col<1,1,1>();
      else return eval_col<1,1,0>();
    } else {
      if (force->newton_pair) return eval_col<1,0,1>();
      else return eval_col<1,0,0>();
    }
  } else {
    if (force->newton_pair) return eval_col<0,0,1>();
    else return eval_col<0,0,0>();
  }
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR> 
void PairYukawaColloidOMP::eval_col()
{

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int i,j,ii,jj,inum,jnum,itype,jtype,tid;
    double xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair,radi,radj;
    double rsq,r2inv,r,rinv,screening,forceyukawa,factor_coul;
    int *ilist,*jlist,*numneigh,**firstneigh;

    ecoul = 0.0;

    double **x = atom->x;
    double **shape = atom->shape;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;
    int nthreads = comm->nthreads;

    double *special_coul = force->special_coul;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  
    // loop over neighbors of my atoms
    int iifrom, iito;
    double **f = loop_setup_thr(atom->f,iifrom,iito,tid,inum,nall,nthreads);

    for (ii = iifrom; ii < iito; ++ii) {

      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      radi = shape[itype][0];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];

	if (j < nall) factor_coul = 1.0;
	else {
	  factor_coul = special_coul[j/nall];
	  j %= nall;
	}

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	jtype = type[j];
	radj = shape[jtype][0];
    
	if (rsq < cutsq[itype][jtype]) {
	  r2inv = 1.0/rsq;
	  r = sqrt(rsq);
	  rinv = 1.0/r;
	  screening = exp(-kappa*(r-(radi+radj)));
	  forceyukawa = a[itype][jtype] * screening;

	  fpair = factor_coul*forceyukawa * rinv;

	  f[i][0] += delx*fpair;
	  f[i][1] += dely*fpair;
	  f[i][2] += delz*fpair;
	  if (NEWTON_PAIR || j < nlocal) {
	    f[j][0] -= delx*fpair;
	    f[j][1] -= dely*fpair;
	    f[j][2] -= delz*fpair;
	  }

	  if (EFLAG) {
	    ecoul = a[itype][jtype]/kappa * screening - offset[itype][jtype];
	    ecoul *= factor_coul;
	  }

	  if (EVFLAG) ev_tally_thr(i,j,nlocal,NEWTON_PAIR,
				   0.0,ecoul,fpair,delx,dely,delz,tid);
	}
      }
    }
    // reduce per thread forces and torques into global force/torque arrays.
    force_reduce_thr(atom->f, nall, nthreads, tid);
  }
  ev_reduce_thr();
  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairYukawaColloidOMP::init_style()
{
  if (!atom->avec->shape_type)
    error->all("Pair yukawa/colloid requires atom attribute shape");
  if (atom->radius_flag)
    error->all("Pair yukawa/colloid cannot be used with "
	       "atom attribute diameter");
  
  // insure all particle shapes are spherical
  // can be point particles or polydisperse

  for (int i = 1; i <= atom->ntypes; i++)
    if ((atom->shape[i][0] != atom->shape[i][1]) || 
	(atom->shape[i][0] != atom->shape[i][2]) ||
	(atom->shape[i][1] != atom->shape[i][2]))
      error->all("Pair yukawa/colloid requires spherical particles");

  int irequest = neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairYukawaColloidOMP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    a[i][j] = mix_energy(a[i][i],a[j][j],1.0,1.0);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  if (offset_flag) {
    double radi = atom->shape[i][0];
    double radj = atom->shape[j][0];
    double screening = exp(-kappa * (cut[i][j] - (radi+radj)));
    offset[i][j] = a[i][j]/kappa * screening;
  } else offset[i][j] = 0.0;

  a[j][i] = a[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairYukawaColloidOMP::single(int i, int j, int itype, int jtype,
				 double rsq,
				 double factor_coul, double factor_lj,
				 double &fforce)
{
  double r2inv,r,rinv,screening,forceyukawa,phi,radi,radj;

  int *type = atom->type;
  radi = atom->shape[itype][0];
  radj = atom->shape[jtype][0];

  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  rinv = 1.0/r;
  screening = exp(-kappa*(r-(radi+radj)));
  forceyukawa = a[itype][jtype] * screening;
  fforce = factor_coul*forceyukawa * rinv;

  phi = a[itype][jtype]/kappa * screening  - offset[itype][jtype]; 
  return factor_coul*phi;
}

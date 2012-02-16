/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "pair_dipole_cut_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairDipoleCutOMP::PairDipoleCutOMP(LAMMPS *lmp) :
  PairDipoleCut(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairDipoleCutOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (evflag) {
      if (eflag) {
	if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
	else eval<1,1,0>(ifrom, ito, thr);
      } else {
	if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
	else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }

    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairDipoleCutOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,qtmp,delx,dely,delz,evdwl,ecoul;
  double rsq,rinv,r2inv,r6inv,r3inv,r5inv,r7inv,fx,fy,fz;
  double forcecoulx,forcecouly,forcecoulz,crossx,crossy,crossz;
  double tixcoul,tiycoul,tizcoul,tjxcoul,tjycoul,tjzcoul;
  double fq,pdotp,pidotr,pjdotr,pre1,pre2,pre3,pre4;
  double forcelj,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  double * const * const torque = thr->get_torque();
  const double * const q = atom->q;
  const double * const * const mu = atom->mu;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  double fxtmp,fytmp,fztmp,t1tmp,t2tmp,t3tmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    qtmp = q[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=t1tmp=t2tmp=t3tmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	r2inv = 1.0/rsq;
	rinv = sqrt(r2inv);

	// atom can have both a charge and dipole
	// i,j = charge-charge, dipole-dipole, dipole-charge, or charge-dipole

	forcecoulx = forcecouly = forcecoulz = 0.0;
	tixcoul = tiycoul = tizcoul = 0.0;
	tjxcoul = tjycoul = tjzcoul = 0.0;
	
	if (rsq < cut_coulsq[itype][jtype]) {

	  if (qtmp != 0.0 && q[j] != 0.0) {
            r3inv = r2inv*rinv;
	    pre1 = qtmp*q[j]*r3inv;

	    forcecoulx += pre1*delx;
	    forcecouly += pre1*dely;
	    forcecoulz += pre1*delz;
	  }

	  if (mu[i][3] > 0.0 && mu[j][3] > 0.0) { 
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
	    r7inv = r5inv*r2inv;

            pdotp = mu[i][0]*mu[j][0] + mu[i][1]*mu[j][1] + mu[i][2]*mu[j][2];
            pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
            pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;

	    pre1 = 3.0*r5inv*pdotp - 15.0*r7inv*pidotr*pjdotr;
	    pre2 = 3.0*r5inv*pjdotr;
	    pre3 = 3.0*r5inv*pidotr;
	    pre4 = -1.0*r3inv;

	    forcecoulx += pre1*delx + pre2*mu[i][0] + pre3*mu[j][0];
	    forcecouly += pre1*dely + pre2*mu[i][1] + pre3*mu[j][1];
	    forcecoulz += pre1*delz + pre2*mu[i][2] + pre3*mu[j][2];
	    
	    crossx = pre4 * (mu[i][1]*mu[j][2] - mu[i][2]*mu[j][1]);
	    crossy = pre4 * (mu[i][2]*mu[j][0] - mu[i][0]*mu[j][2]);
	    crossz = pre4 * (mu[i][0]*mu[j][1] - mu[i][1]*mu[j][0]);

	    tixcoul += crossx + pre2 * (mu[i][1]*delz - mu[i][2]*dely);
	    tiycoul += crossy + pre2 * (mu[i][2]*delx - mu[i][0]*delz);
	    tizcoul += crossz + pre2 * (mu[i][0]*dely - mu[i][1]*delx);
	    tjxcoul += -crossx + pre3 * (mu[j][1]*delz - mu[j][2]*dely);
	    tjycoul += -crossy + pre3 * (mu[j][2]*delx - mu[j][0]*delz);
	    tjzcoul += -crossz + pre3 * (mu[j][0]*dely - mu[j][1]*delx);
	  }

	  if (mu[i][3] > 0.0 && q[j] != 0.0) { 
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
	    pre1 = 3.0*q[j]*r5inv * pidotr;
	    pre2 = q[j]*r3inv;

	    forcecoulx += pre2*mu[i][0] - pre1*delx;
            forcecouly += pre2*mu[i][1] - pre1*dely;
            forcecoulz += pre2*mu[i][2] - pre1*delz;
	    tixcoul += pre2 * (mu[i][1]*delz - mu[i][2]*dely);
	    tiycoul += pre2 * (mu[i][2]*delx - mu[i][0]*delz);
	    tizcoul += pre2 * (mu[i][0]*dely - mu[i][1]*delx);
	  }

	  if (mu[j][3] > 0.0 && qtmp != 0.0) { 
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;
	    pre1 = 3.0*qtmp*r5inv * pjdotr;
	    pre2 = qtmp*r3inv;

	    forcecoulx += pre1*delx - pre2*mu[j][0];
            forcecouly += pre1*dely - pre2*mu[j][1];
            forcecoulz += pre1*delz - pre2*mu[j][2];
	    tjxcoul += -pre2 * (mu[j][1]*delz - mu[j][2]*dely);
	    tjycoul += -pre2 * (mu[j][2]*delx - mu[j][0]*delz);
	    tjzcoul += -pre2 * (mu[j][0]*dely - mu[j][1]*delx);
	  }
	}

	// LJ interaction

	if (rsq < cut_ljsq[itype][jtype]) {
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  forcelj *= factor_lj * r2inv;
	} else forcelj = 0.0;
	  
	// total force

	fq = factor_coul*qqrd2e;
	fx = fq*forcecoulx + delx*forcelj;
	fy = fq*forcecouly + dely*forcelj;
	fz = fq*forcecoulz + delz*forcelj;
	
	// force & torque accumulation

	fxtmp += fx;
	fytmp += fy;
	fztmp += fz;
	t1tmp += fq*tixcoul;
	t2tmp += fq*tiycoul;
	t3tmp += fq*tizcoul;

	if (NEWTON_PAIR || j < nlocal) {
	  f[j][0] -= fx;
	  f[j][1] -= fy;
	  f[j][2] -= fz;
	  torque[j][0] += fq*tjxcoul;
	  torque[j][1] += fq*tjycoul;
	  torque[j][2] += fq*tjzcoul;
	}

	if (EFLAG) {
	  if (rsq < cut_coulsq[itype][jtype]) {
	    ecoul = qtmp*q[j]*rinv;
	    if (mu[i][3] > 0.0 && mu[j][3] > 0.0)
	      ecoul += r3inv*pdotp - 3.0*r5inv*pidotr*pjdotr;
	    if (mu[i][3] > 0.0 && q[j] != 0.0) 
	      ecoul += -q[j]*r3inv*pidotr;
	    if (mu[j][3] > 0.0 && qtmp != 0.0)
	      ecoul += qtmp*r3inv*pjdotr;
	    ecoul *= factor_coul*qqrd2e;
	  } else ecoul = 0.0;

	  if (rsq < cut_ljsq[itype][jtype]) {
	    evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
	      offset[itype][jtype];
	    evdwl *= factor_lj;
	  } else evdwl = 0.0;
	}

	if (EVFLAG) ev_tally_xyz_thr(this,i,j,nlocal,NEWTON_PAIR,
				     evdwl,ecoul,fx,fy,fz,delx,dely,delz,thr);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
    torque[i][0] += t1tmp;
    torque[i][1] += t2tmp;
    torque[i][2] += t3tmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairDipoleCutOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairDipoleCut::memory_usage();

  return bytes;
}

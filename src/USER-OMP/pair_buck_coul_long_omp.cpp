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
#include "pair_buck_coul_long_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairBuckCoulLongOMP::PairBuckCoulLongOMP(LAMMPS *lmp) :
  PairBuckCoulLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulLongOMP::compute(int eflag, int vflag)
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

    // reduce per thread forces into global force array.
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairBuckCoulLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,r,rexp,forcecoul,forcebuck,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;
  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	r2inv = 1.0/rsq;
	r = sqrt(rsq);

	if (rsq < cut_coulsq) {
	  grij = g_ewald * r;
	  expm2 = exp(-grij*grij);
	  t = 1.0 / (1.0 + EWALD_P*grij);
	  erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
	  prefactor = qqrd2e * qtmp*q[j]/r;
	  forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	  if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
	} else forcecoul = 0.0;

	if (rsq < cut_ljsq[itype][jtype]) {
	  r6inv = r2inv*r2inv*r2inv;
	  rexp = exp(-r*rhoinv[itype][jtype]);
	  forcebuck = buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv;
	} else forcebuck = 0.0;
	
	fpair = (forcecoul + factor_lj*forcebuck)*r2inv;

	fxtmp += delx*fpair;
	fytmp += dely*fpair;
	fztmp += delz*fpair;
	if (NEWTON_PAIR || j < nlocal) {
	  f[j][0] -= delx*fpair;
	  f[j][1] -= dely*fpair;
	  f[j][2] -= delz*fpair;
	}

	if (EFLAG) {
	  if (rsq < cut_coulsq) {
	    ecoul = prefactor*erfc;
	    if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
	  } else ecoul = 0.0;
	  if (rsq < cut_ljsq[itype][jtype]) {
	    evdwl = a[itype][jtype]*rexp - c[itype][jtype]*r6inv -
	      offset[itype][jtype];
	    evdwl *= factor_lj;
	  } else evdwl = 0.0;
	}

	if (EVFLAG) ev_tally_thr(this, i,j,nlocal,NEWTON_PAIR,
				 evdwl,ecoul,fpair,delx,dely,delz,thr);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairBuckCoulLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBuckCoulLong::memory_usage();

  return bytes;
}

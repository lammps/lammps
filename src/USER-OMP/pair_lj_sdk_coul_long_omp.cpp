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
#include "pair_lj_sdk_coul_long_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "lj_sdk_common.h"

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace LJSDKParms;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJSDKCoulLongOMP::PairLJSDKCoulLongOMP(LAMMPS *lmp) :
  PairLJSDKCoulLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJSDKCoulLongOMP::compute(int eflag, int vflag)
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
	if (force->newton_pair) eval_thr<1,1,1>(ifrom, ito, thr);
	else eval_thr<1,1,0>(ifrom, ito, thr);
      } else {
	if (force->newton_pair) eval_thr<1,0,1>(ifrom, ito, thr);
	else eval_thr<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval_thr<0,0,1>(ifrom, ito, thr);
      else eval_thr<0,0,0>(ifrom, ito, thr);
    }

    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJSDKCoulLongOMP::eval_thr(int iifrom, int iito, ThrData * const thr)
{
  int i,ii,j,jj,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double r,rsq,r2inv,forcecoul,forcelj,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  double fxtmp,fytmp,fztmp;

  const int * const ilist = list->ilist;
  const int * const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp=fytmp=fztmp=0.0;

    const int itype = type[i];
    const int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      forcecoul = forcelj = evdwl = ecoul = 0.0;

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
	const int ljt = lj_type[itype][jtype];

	if (rsq < cut_coulsq) {
	  if (!ncoultablebits || rsq <= tabinnersq) {
	    r = sqrt(rsq);
	    grij = g_ewald * r;
	    expm2 = exp(-grij*grij);
	    t = 1.0 / (1.0 + EWALD_P*grij);
	    erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
	    prefactor = qqrd2e * qtmp*q[j]/r;
	    forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	    if (EFLAG)
	      ecoul = prefactor*erfc;
	    if (factor_coul < 1.0) {
	      forcecoul -= (1.0-factor_coul)*prefactor;
	      if (EFLAG)
		ecoul -= (1.0-factor_coul)*prefactor;
	    }
	  } else {
	    union_int_float_t rsq_lookup;
	    rsq_lookup.f = rsq;
	    itable = rsq_lookup.i & ncoulmask;
	    itable >>= ncoulshiftbits;
	    fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
	    table = ftable[itable] + fraction*dftable[itable];
	    forcecoul = qtmp*q[j] * table;
	    if (EFLAG)
	      ecoul = qtmp*q[j] * table;
	    if (factor_coul < 1.0) {
	      table = ctable[itable] + fraction*dctable[itable];
	      prefactor = qtmp*q[j] * table;
	      forcecoul -= (1.0-factor_coul)*prefactor;
	      if (EFLAG)
		ecoul -= (1.0-factor_coul)*prefactor;
	    }
	  }
	}

	if (rsq < cut_ljsq[itype][jtype]) {

	  if (ljt == LJ12_4) {
	    const double r4inv=r2inv*r2inv;
	    forcelj = r4inv*(lj1[itype][jtype]*r4inv*r4inv
			     - lj2[itype][jtype]);

	    if (EFLAG)
	      evdwl = r4inv*(lj3[itype][jtype]*r4inv*r4inv
			     - lj4[itype][jtype]) - offset[itype][jtype];

	  } else if (ljt == LJ9_6) {
	    const double r3inv = r2inv*sqrt(r2inv);
	    const double r6inv = r3inv*r3inv;
	    forcelj = r6inv*(lj1[itype][jtype]*r3inv
			     - lj2[itype][jtype]);
	    if (EFLAG)
	      evdwl = r6inv*(lj3[itype][jtype]*r3inv
			     - lj4[itype][jtype]) - offset[itype][jtype];

	  } else if (ljt == LJ12_6) {
	    const double r6inv = r2inv*r2inv*r2inv;
	    forcelj = r6inv*(lj1[itype][jtype]*r6inv
			     - lj2[itype][jtype]);
	    if (EFLAG)
	      evdwl = r6inv*(lj3[itype][jtype]*r6inv
			     - lj4[itype][jtype]) - offset[itype][jtype];
	  }
	  forcelj *= factor_lj;
	  if (EFLAG)
	    evdwl *= factor_lj;
	}

	fpair = (forcecoul + forcelj) * r2inv;

	fxtmp += delx*fpair;
	fytmp += dely*fpair;
	fztmp += delz*fpair;
	if (NEWTON_PAIR || j < nlocal) {
	  f[j][0] -= delx*fpair;
	  f[j][1] -= dely*fpair;
	  f[j][2] -= delz*fpair;
	}

	if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
				 evdwl,ecoul,fpair,delx,dely,delz,thr);

      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJSDKCoulLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJSDKCoulLong::memory_usage();

  return bytes;
}

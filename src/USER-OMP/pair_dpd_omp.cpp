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
#include "pair_dpd_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "random_mars.h"

using namespace LAMMPS_NS;

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairDPDOMP::PairDPDOMP(LAMMPS *lmp) :
  PairDPD(lmp), ThrOMP(lmp, PAIR)
{
  respa_enable = 0;
  random_thr = NULL;
}

/* ---------------------------------------------------------------------- */

PairDPDOMP::~PairDPDOMP() 
{
  if (random_thr) {
    for (int i=1; i < comm->nthreads; ++i)
      delete random_thr[i];

    delete[] random_thr;
    random_thr = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void PairDPDOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(this);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  if (!random_thr)
    random_thr = new RanMars*[nthreads];
  
  random_thr[0] = random;

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int ifrom, ito, tid;
    double **f;

    f = loop_setup_thr(atom->f, ifrom, ito, tid, inum, nall, nthreads);

    if (random_thr && tid > 0)
      random_thr[tid] = new RanMars(Pair::lmp, seed + comm->me 
				    + comm->nprocs*tid);

    if (evflag) {
      if (eflag) {
	if (force->newton_pair) eval<1,1,1>(f, ifrom, ito, tid);
	else eval<1,1,0>(f, ifrom, ito, tid);
      } else {
	if (force->newton_pair) eval<1,0,1>(f, ifrom, ito, tid);
	else eval<1,0,0>(f, ifrom, ito, tid);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(f, ifrom, ito, tid);
      else eval<0,0,0>(f, ifrom, ito, tid);
    }

    // reduce per thread forces into global force array.
    data_reduce_thr(&(atom->f[0][0]), nall, nthreads, 3, tid);
  } // end of omp parallel region

  // reduce per thread energy and virial, if requested.
  if (evflag) ev_reduce_thr(this);
  if (vflag_fdotr) virial_fdotr_compute();
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairDPDOMP::eval(double **f, int iifrom, int iito, int tid)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,randnum,factor_dpd;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  double dtinvsqrt = 1.0/sqrt(update->dt);
  double fxtmp,fytmp,fztmp;
  RanMars &rng = *random_thr[tid];

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	r = sqrt(rsq);
	if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
	rinv = 1.0/r;
	delvx = vxtmp - v[j][0];
	delvy = vytmp - v[j][1];
	delvz = vztmp - v[j][2];
	dot = delx*delvx + dely*delvy + delz*delvz;
	wd = 1.0 - r/cut[itype][jtype];
	randnum = rng.gaussian();

	// conservative force = a0 * wd
	// drag force = -gamma * wd^2 * (delx dot delv) / r
	// random force = sigma * wd * rnd * dtinvsqrt;

	fpair = a0[itype][jtype]*wd;
	fpair -= gamma[itype][jtype]*wd*wd*dot*rinv;
	fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
	fpair *= factor_dpd*rinv;	

	fxtmp += delx*fpair;
	fytmp += dely*fpair;
	fztmp += delz*fpair;
	if (NEWTON_PAIR || j < nlocal) {
	  f[j][0] -= delx*fpair;
	  f[j][1] -= dely*fpair;
	  f[j][2] -= delz*fpair;
	}

	if (EFLAG) {
          // unshifted eng of conservative term:
	  // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
	  // eng shifted to 0.0 at cutoff
	  evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
	  evdwl *= factor_dpd;
	}

	if (EVFLAG) ev_tally_thr(this, i,j,nlocal,NEWTON_PAIR,
				 evdwl,0.0,fpair,delx,dely,delz,tid);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairDPDOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairDPD::memory_usage();
  bytes += comm->nthreads * sizeof(RanMars*);
  bytes += comm->nthreads * sizeof(RanMars);

  return bytes;
}

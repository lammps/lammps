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
#include "pair_hbond_dreiding_morse_omp.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "math_const.h"

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

PairHbondDreidingMorseOMP::PairHbondDreidingMorseOMP(LAMMPS *lmp) :
  PairHbondDreidingMorse(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  hbcount_thr = hbeng_thr = NULL;
}

/* ---------------------------------------------------------------------- */

PairHbondDreidingMorseOMP::~PairHbondDreidingMorseOMP()
{
  if (hbcount_thr) {
    delete[] hbcount_thr;
    delete[] hbeng_thr;
  }
}

/* ---------------------------------------------------------------------- */

void PairHbondDreidingMorseOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  if (!hbcount_thr) {
    hbcount_thr = new double[nthreads];
    hbeng_thr = new double[nthreads];
  }

  for (int i=0; i < nthreads; ++i) {
    hbcount_thr[i] = 0.0;
    hbeng_thr[i] = 0.0;
  }

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

  // reduce per thread hbond data
  if (eflag_global) {
    pvector[0] = 0.0;
    pvector[1] = 0.0;
    for (int i=0; i < nthreads; ++i) {
      pvector[0] += hbcount_thr[i];
      pvector[1] += hbeng_thr[i];
    }
  }
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairHbondDreidingMorseOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,k,m,ii,jj,kk,jnum,itype,jtype,ktype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,rsq1,rsq2,r1,r2;
  double factor_hb,force_angle,force_kernel,evdwl;
  double c,s,a,b,ac,a11,a12,a22,vx1,vx2,vy1,vy2,vz1,vz2;
  double fi[3],fj[3],delr1[3],delr2[3];
  double r,dr,dexp,eng_morse,switch1,switch2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  Param *pm;

  evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const type = atom->type;
  const double * const special_lj = force->special_lj;
  const int * const * const nspecial = atom->nspecial;
  const int * const * const special = atom->special;
  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  
  // ii = loop over donors
  // jj = loop over acceptors
  // kk = loop over hydrogens bonded to donor

  int hbcount = 0;
  double hbeng = 0.0;

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    itype = type[i];
    if (!donor[itype]) continue;

    const int * const klist = special[i];
    const int knum = nspecial[i][0];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_hb = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      jtype = type[j];
      if (!acceptor[jtype]) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      for (kk = 0; kk < knum; kk++) {
	k = atom->map(klist[kk]);
	if (k < 0) continue;
	ktype = type[k];
	m = type2param[itype][jtype][ktype];
	if (m < 0) continue;
	pm = &params[m];

	if (rsq < pm->cut_outersq) {
	  delr1[0] = xtmp - x[k][0];
	  delr1[1] = ytmp - x[k][1];
	  delr1[2] = ztmp - x[k][2];
	  domain->minimum_image(delr1);
	  rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
	  r1 = sqrt(rsq1);
	  
	  delr2[0] = x[j][0] - x[k][0];
	  delr2[1] = x[j][1] - x[k][1];
	  delr2[2] = x[j][2] - x[k][2];
	  domain->minimum_image(delr2);
	  rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
	  r2 = sqrt(rsq2);
	  
	  // angle (cos and sin)
	  
	  c = delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2];
	  c /= r1*r2;
	  if (c > 1.0) c = 1.0;
	  if (c < -1.0) c = -1.0;
	  ac = acos(c);

	  if (ac > pm->cut_angle && ac < (2.0*MY_PI - pm->cut_angle)) {
	    s = sqrt(1.0 - c*c);
	    if (s < SMALL) s = SMALL;

	    // Morse-specific kernel

	    r = sqrt(rsq);
	    dr = r - pm->r0;
	    dexp = exp(-pm->alpha * dr);
	    force_kernel = pm->morse1*(dexp*dexp - dexp)/r * pow(c,pm->ap);
	    force_angle = pm->ap * eng_morse * pow(c,pm->ap-1)*s;

	    eng_morse = pm->d0 * (dexp*dexp - 2.0*dexp);
	    if (rsq > pm->cut_innersq) {
	      switch1 = (pm->cut_outersq-rsq) * (pm->cut_outersq-rsq) *
			(pm->cut_outersq + 2.0*rsq - 3.0*pm->cut_innersq) /
			pm->denom_vdw;
	      switch2 = 12.0*rsq * (pm->cut_outersq-rsq) *
			(rsq-pm->cut_innersq) / pm->denom_vdw;
	      force_kernel = force_kernel*switch1 + eng_morse*switch2;
	      eng_morse *= switch1;
	    }

	    if (EFLAG) {
	      evdwl = eng_morse * pow(c,params[m].ap);
	      evdwl *= factor_hb;
	    }

	    a = factor_hb*force_angle/s;
	    b = factor_hb*force_kernel;
	    
	    a11 = a*c / rsq1;
	    a12 = -a / (r1*r2);
	    a22 = a*c / rsq2;
	    
	    vx1 = a11*delr1[0] + a12*delr2[0];
	    vx2 = a22*delr2[0] + a12*delr1[0];
	    vy1 = a11*delr1[1] + a12*delr2[1];
	    vy2 = a22*delr2[1] + a12*delr1[1];
	    vz1 = a11*delr1[2] + a12*delr2[2];
	    vz2 = a22*delr2[2] + a12*delr1[2];
	    
	    fi[0] = vx1 + b*delx;
	    fi[1] = vy1 + b*dely;
	    fi[2] = vz1 + b*delz;
	    fj[0] = vx2 - b*delx;
	    fj[1] = vy2 - b*dely;
	    fj[2] = vz2 - b*delz;

	    fxtmp += fi[0];
	    fytmp += fi[1];
	    fztmp += fi[2];

	    f[j][0] += fj[0];
	    f[j][1] += fj[1];
	    f[j][2] += fj[2];
	    
	    f[k][0] -= vx1 + vx2;
	    f[k][1] -= vy1 + vy2;
	    f[k][2] -= vz1 + vz2;

	    // KIJ instead of IJK b/c delr1/delr2 are both with respect to k

	    if (EVFLAG) ev_tally3_thr(this,k,i,j,evdwl,0.0,fi,fj,delr1,delr2,thr);
	    if (EFLAG) {
	      hbcount++;
	      hbeng += evdwl;
	    }
	  }
	}
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
  const int tid = thr->get_tid();
  hbcount_thr[tid] = static_cast<double>(hbcount);
  hbeng_thr[tid] = hbeng;
}

/* ---------------------------------------------------------------------- */

double PairHbondDreidingMorseOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += comm->nthreads * 2 * sizeof(double);
  bytes += PairHbondDreidingMorse::memory_usage();

  return bytes;
}

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
#include "pair_buck_coul_omp.h"
#include "atom.h"
#include "comm.h"
#include "math_vector.h"
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

PairBuckCoulOMP::PairBuckCoulOMP(LAMMPS *lmp) :
  PairBuckCoul(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulOMP::compute(int eflag, int vflag)
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

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairBuckCoulOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  double evdwl,ecoul,fpair;
  evdwl = ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;

  const double *x0 = x[0];
  double *f0 = f[0], *fi = f0;

  int *ilist = list->ilist;

  // loop over neighbors of my atoms

  int i, ii, j, order1 = ewald_order&(1<<1), order6 = ewald_order&(1<<6);
  int *jneigh, *jneighn, typei, typej, ni;
  double qi, qri, *cutsqi, *cut_bucksqi,
	 *buck1i, *buck2i, *buckai, *buckci, *rhoinvi, *offseti;
  double r, rsq, r2inv, force_coul, force_buck;
  double g2 = g_ewald*g_ewald, g6 = g2*g2*g2, g8 = g6*g2;
  vector xi, d;

  for (ii = iifrom; ii < iito; ++ii) {			// loop over my atoms
    i = ilist[ii]; fi = f0+3*i;
    if (order1) qri = (qi = q[i])*qqrd2e;		// initialize constants
    offseti = offset[typei = type[i]];
    buck1i = buck1[typei]; buck2i = buck2[typei];
    buckai = buck_a[typei]; buckci = buck_c[typei], rhoinvi = rhoinv[typei];
    cutsqi = cutsq[typei]; cut_bucksqi = cut_bucksq[typei];
    memcpy(xi, x0+(i+(i<<1)), sizeof(vector));
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {			// loop over neighbors
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;
      
      { const register double *xj = x0+(j+(j<<1));
	d[0] = xi[0] - xj[0];				// pair vector
	d[1] = xi[1] - xj[1];
	d[2] = xi[2] - xj[2]; }

      if ((rsq = vec_dot(d, d)) >= cutsqi[typej = type[j]]) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      if (order1 && (rsq < cut_coulsq)) {		// coulombic
	if (!ncoultablebits || rsq <= tabinnersq) {	// series real space
	  register double x = g_ewald*r;
	  register double s = qri*q[j], t = 1.0/(1.0+EWALD_P*x);
	  if (ni == 0) {
	    s *= g_ewald*exp(-x*x);
	    force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s;
	    if (EFLAG) ecoul = t;
	  } else {					// special case
	    register double f = s*(1.0-special_coul[ni])/r;
	    s *= g_ewald*exp(-x*x);
	    force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-f;
	    if (EFLAG) ecoul = t-f;
	  }					// table real space
	} else {
	  register union_int_float_t t;
	  t.f = rsq;
	  register const int k = (t.i & ncoulmask) >> ncoulshiftbits;
	  register double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
	  if (ni == 0) {
	    force_coul = qiqj*(ftable[k]+f*dftable[k]);
	    if (EFLAG) ecoul = qiqj*(etable[k]+f*detable[k]);
	  }
	  else {					// special case
	    t.f = (1.0-special_coul[ni])*(ctable[k]+f*dctable[k]);
	    force_coul = qiqj*(ftable[k]+f*dftable[k]-t.f);
	    if (EFLAG) ecoul = qiqj*(etable[k]+f*detable[k]-t.f);
	  }
	}
      } else force_coul = ecoul = 0.0;

      if (rsq < cut_bucksqi[typej]) {			// buckingham
	register double rn = r2inv*r2inv*r2inv, 
			expr = exp(-r*rhoinvi[typej]);
	if (order6) {					// long-range
	  register double x2 = g2*rsq, a2 = 1.0/x2;
	  x2 = a2*exp(-x2)*buckci[typej];
	  if (ni == 0) {
	    force_buck =
	      r*expr*buck1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
	    if (EFLAG) evdwl = expr*buckai[typej]-g6*((a2+1.0)*a2+0.5)*x2;
	  } else {					// special case
	    register double f = special_lj[ni], t = rn*(1.0-f);
	    force_buck = f*r*expr*buck1i[typej]-
	      g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*buck2i[typej];
	    if (EFLAG) evdwl = f*expr*buckai[typej] - 
			 g6*((a2+1.0)*a2+0.5)*x2+t*buckci[typej];
	  }
	} else {						// cut
	  if (ni == 0) {
	    force_buck = r*expr*buck1i[typej]-rn*buck2i[typej];
	    if (EFLAG) evdwl = expr*buckai[typej] - 
			 rn*buckci[typej]-offseti[typej];
	  } else {					// special case
	    register double f = special_lj[ni];
	    force_buck = f*(r*expr*buck1i[typej]-rn*buck2i[typej]);
	    if (EFLAG) 
	      evdwl = f*(expr*buckai[typej]-rn*buckci[typej]-offseti[typej]);
	  }
	}
      } else force_buck = evdwl = 0.0;

      fpair = (force_coul+force_buck)*r2inv;

      if (NEWTON_PAIR || j < nlocal) {
	register double *fj = f0+(j+(j<<1)), f;
	fi[0] += f = d[0]*fpair; fj[0] -= f;
	fi[1] += f = d[1]*fpair; fj[1] -= f;
	fi[2] += f = d[2]*fpair; fj[2] -= f;
      } else {
	fi[0] += d[0]*fpair;
	fi[1] += d[1]*fpair;
	fi[2] += d[2]*fpair;
      }

      if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
			       evdwl,ecoul,fpair,d[0],d[1],d[2],thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairBuckCoulOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBuckCoul::memory_usage();

  return bytes;
}

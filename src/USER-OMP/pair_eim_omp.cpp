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
#include "string.h"

#include "pair_eim_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEIMOMP::PairEIMOMP(LAMMPS *lmp) :
  PairEIM(lmp), ThrOMP(lmp, PAIR)
{
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairEIMOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(this);
  } else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho,nthreads*nmax,"pair:rho");
    memory->create(fp,nthreads*nmax,"pair:fp");
  }

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int ifrom, ito, tid;
    double **f, *rho_t, *fp_t;

    f = loop_setup_thr(atom->f, ifrom, ito, tid, inum, nall, nthreads);
    if (force->newton_pair) {
      rho_t = rho + tid*nall;
      fp_t = fp + tid*nall;
    } else {
      rho_t = rho + tid*atom->nlocal;
      fp_t = fp + tid*atom->nlocal;
    }
    
    if (evflag) {
      if (eflag) {
	if (force->newton_pair) eval<1,1,1>(f, rho_t, fp_t, ifrom, ito, tid);
	else eval<1,1,0>(f, rho_t, fp_t, ifrom, ito, tid);
      } else {
	if (force->newton_pair) eval<1,0,1>(f, rho_t, fp_t, ifrom, ito, tid);
	else eval<1,0,0>(f, rho_t, fp_t, ifrom, ito, tid);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(f, rho_t, fp_t, ifrom, ito, tid);
      else eval<0,0,0>(f, rho_t, fp_t, ifrom, ito, tid);
    }

    // reduce per thread forces into global force array.
    data_reduce_thr(&(atom->f[0][0]), nall, nthreads, 3, tid);
  } // end of omp parallel region

  // reduce per thread energy and virial, if requested.
  if (evflag) ev_reduce_thr(this);
  if (vflag_fdotr) virial_fdotr_compute();
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairEIMOMP::eval(double **f, double *rho_t, double *fp_t,
		      int iifrom, int iito, int tid)
{
  int i,j,ii,jj,m,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,phip,phi,coul,coulp,recip,psip;
  double *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density and fp

  if (NEWTON_PAIR) {
    memset(rho_t, 0, nall*sizeof(double));
    memset(fp_t, 0, nall*sizeof(double));
  } else {
    memset(rho_t, 0, nlocal*sizeof(double));
    memset(fp_t, 0, nlocal*sizeof(double));
  }

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq[itype][jtype]) {
	p = sqrt(rsq)*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);
	coeff = Fij_spline[type2Fij[itype][jtype]][m];
	rho_t[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
	if (NEWTON_PAIR || j < nlocal) {
	  coeff = Fij_spline[type2Fij[jtype][itype]][m];
	  rho_t[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
	}
      }
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // communicate and sum densities
  if (NEWTON_PAIR) {
    // reduce per thread density
    data_reduce_thr(&(rho[0]), nall, comm->nthreads, 1, tid);

    // wait until reduction is complete
    sync_threads();

#if defined(_OPENMP)
#pragma omp master
#endif
    { 
      rhofp = 1;
      comm->reverse_comm_pair(this); 
    }

  } else {
    data_reduce_thr(&(rho[0]), nlocal, comm->nthreads, 1, tid);

    // wait until reduction is complete
    sync_threads();
  }

#if defined(_OPENMP)
#pragma omp master
#endif
  { 
    rhofp = 1;
    comm->forward_comm_pair(this); 
  }

  // wait until master is finished communicating
  sync_threads();
 
  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
 
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
 
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
 
      if (rsq < cutforcesq[itype][jtype]) {
        p = sqrt(rsq)*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = Gij_spline[type2Gij[itype][jtype]][m];
        fp_t[i] += rho[j]*(((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6]);
        if (NEWTON_PAIR || j < nlocal) {
          fp_t[j] += rho[i]*(((coeff[3]*p + coeff[4])*p + coeff[5])*p + 
			   coeff[6]);
        }
      }
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // communicate and sum modified densities
  if (NEWTON_PAIR) {
    // reduce per thread density
    data_reduce_thr(&(fp[0]), nall, comm->nthreads, 1, tid);

    // wait until reduction is complete
    sync_threads();

#if defined(_OPENMP)
#pragma omp master
#endif
    { 
      rhofp = 2;
      comm->reverse_comm_pair(this); 
    }

  } else {
    data_reduce_thr(&(fp[0]), nlocal, comm->nthreads, 1, tid);

    // wait until reduction is complete
    sync_threads();
  }

#if defined(_OPENMP)
#pragma omp master
#endif
  { 
    rhofp = 2;
    comm->forward_comm_pair(this); 
  }

  // wait until master is finished communicating
  sync_threads();

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    itype = type[i];
    if (EFLAG) {
      phi = 0.5*rho[i]*fp[i];
      if (eflag_global) eng_vdwl_thr[tid] += phi;
      if (eflag_atom) eatom_thr[tid][i] += phi;
    }
  }

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    fxtmp = fytmp = fztmp = 0.0;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq[itype][jtype]) {
	r = sqrt(rsq);
	p = r*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'

        coeff = Fij_spline[type2Fij[jtype][itype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = Fij_spline[type2Fij[itype][jtype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = phiij_spline[type2phiij[itype][jtype]][m];
        phip = (coeff[0]*p + coeff[1])*p + coeff[2];
        phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = Gij_spline[type2Gij[itype][jtype]][m];
        coul = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coulp = (coeff[0]*p + coeff[1])*p + coeff[2];
        psip = phip + (rho[i]*rho[j]-q0[itype]*q0[jtype])*coulp +
               fp[i]*rhojp + fp[j]*rhoip;
        recip = 1.0/r;
        fpair = -psip*recip;
	fxtmp += delx*fpair;
	fytmp += dely*fpair;
	fztmp += delz*fpair;
	if (NEWTON_PAIR || j < nlocal) {
	  f[j][0] -= delx*fpair;
	  f[j][1] -= dely*fpair;
	  f[j][2] -= delz*fpair;
	}

	if (EFLAG) evdwl = phi-q0[itype]*q0[jtype]*coul;
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

double PairEIMOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairEIM::memory_usage();

  return bytes;
}

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

#include "pair_adp_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairADPOMP::PairADPOMP(LAMMPS *lmp) :
  PairADP(lmp), ThrOMP(lmp, PAIR)
{
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairADPOMP::compute(int eflag, int vflag)
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
    memory->destroy(mu);
    memory->destroy(lambda);
    nmax = atom->nmax;
    memory->create(rho,nthreads*nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
    memory->create(mu,nthreads*nmax,3,"pair:mu");
    memory->create(lambda,nthreads*nmax,6,"pair:lambda");
  }

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int ifrom, ito, tid;
    double **f, *rho_t, **mu_t, **lambda_t;

    f = loop_setup_thr(atom->f, ifrom, ito, tid, inum, nall, nthreads);
    if (force->newton_pair) {
      rho_t = rho + tid*nall;
      mu_t = mu + tid*nall;
      lambda_t = lambda + tid*nall;
    } else {
      rho_t = rho + tid*atom->nlocal;
      mu_t = mu + tid*atom->nlocal;
      lambda_t = lambda + tid*atom->nlocal;
    }

    if (evflag) {
      if (eflag) {
	if (force->newton_pair) eval<1,1,1>(f, rho_t, mu_t, lambda_t, ifrom, ito, tid);
	else eval<1,1,0>(f, rho_t, mu_t, lambda_t, ifrom, ito, tid);
      } else {
	if (force->newton_pair) eval<1,0,1>(f, rho_t, mu_t, lambda_t, ifrom, ito, tid);
	else eval<1,0,0>(f, rho_t, mu_t, lambda_t, ifrom, ito, tid);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(f, rho_t, mu_t, lambda_t, ifrom, ito, tid);
      else eval<0,0,0>(f, rho_t, mu_t, lambda_t, ifrom, ito, tid);
    }

    // reduce per thread forces into global force array.
    data_reduce_thr(&(atom->f[0][0]), nall, nthreads, 3, tid);
  } // end of omp parallel region

  // reduce per thread energy and virial, if requested.
  if (evflag) ev_reduce_thr(this);
  if (vflag_fdotr) virial_fdotr_compute();
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairADPOMP::eval(double **f, double *rho_t, double **mu_t, 
		      double **lambda_t, int iifrom, int iito, int tid)
{
  int i,j,ii,jj,m,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double u2,u2p,w2,w2p,nu;
  double *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double delmux,delmuy,delmuz,trdelmu,tradellam;
  double adpx,adpy,adpz,fx,fy,fz;
  double sumlamxx,sumlamyy,sumlamzz,sumlamyz,sumlamxz,sumlamxy;

  evdwl = 0.0;

  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density 

  if (NEWTON_PAIR) {
    memset(rho_t, 0, nall*sizeof(double));
    memset(&(mu_t[0][0]), 0, 3*nall*sizeof(double));
    memset(&(lambda_t[0][0]), 0, 6*nall*sizeof(double));
  } else {
    memset(rho_t, 0, nlocal*sizeof(double));
    memset(&(mu_t[0][0]), 0, 3*nlocal*sizeof(double));
    memset(&(lambda_t[0][0]), 0, 6*nlocal*sizeof(double));
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

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
	jtype = type[j];
	p = sqrt(rsq)*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);
	coeff = rhor_spline[type2rhor[jtype][itype]][m];
	rho_t[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = u2r_spline[type2u2r[jtype][itype]][m];
	u2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        mu_t[i][0] += u2*delx;
        mu_t[i][1] += u2*dely;
        mu_t[i][2] += u2*delz;
        coeff = w2r_spline[type2w2r[jtype][itype]][m];
	w2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        lambda_t[i][0] += w2*delx*delx;
        lambda_t[i][1] += w2*dely*dely;
        lambda_t[i][2] += w2*delz*delz;
        lambda_t[i][3] += w2*dely*delz;
        lambda_t[i][4] += w2*delx*delz;
        lambda_t[i][5] += w2*delx*dely;

	if (NEWTON_PAIR || j < nlocal) {
          // verify sign difference for mu and lambda
	  coeff = rhor_spline[type2rhor[itype][jtype]][m];
	  rho_t[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          coeff = u2r_spline[type2u2r[itype][jtype]][m];
          u2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          mu_t[j][0] -= u2*delx;
          mu_t[j][1] -= u2*dely;
          mu_t[j][2] -= u2*delz;
          coeff = w2r_spline[type2w2r[itype][jtype]][m];
	  w2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          lambda_t[j][0] += w2*delx*delx;
          lambda_t[j][1] += w2*dely*dely;
          lambda_t[j][2] += w2*delz*delz;
          lambda_t[j][3] += w2*dely*delz;
          lambda_t[j][4] += w2*delx*delz;
          lambda_t[j][5] += w2*delx*dely;
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
    data_reduce_thr(&(mu[0][0]), nall, comm->nthreads, 3, tid);
    data_reduce_thr(&(lambda[0][0]), nall, comm->nthreads, 6, tid);

    // wait until reduction is complete
    sync_threads();

#if defined(_OPENMP)
#pragma omp master
#endif
    { comm->reverse_comm_pair(this); }

    // wait until master thread is done with communication
    sync_threads();
  
  } else {
    // reduce per thread density
    data_reduce_thr(&(rho[0]), nlocal, comm->nthreads, 1, tid);
    data_reduce_thr(&(mu[0][0]), nlocal, comm->nthreads, 3, tid);
    data_reduce_thr(&(lambda[0][0]), nlocal, comm->nthreads, 6, tid);

    // wait until reduction is complete
    sync_threads();
  }
  
  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    p = rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    if (EFLAG) {
      phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      phi += 0.5*(mu[i][0]*mu[i][0]+mu[i][1]*mu[i][1]+mu[i][2]*mu[i][2]);
      phi += 0.5*(lambda[i][0]*lambda[i][0]+lambda[i][1]*
		  lambda[i][1]+lambda[i][2]*lambda[i][2]);
      phi += 1.0*(lambda[i][3]*lambda[i][3]+lambda[i][4]*
		  lambda[i][4]+lambda[i][5]*lambda[i][5]);
      phi -= 1.0/6.0*(lambda[i][0]+lambda[i][1]+lambda[i][2])*
	(lambda[i][0]+lambda[i][1]+lambda[i][2]);
      if (eflag_global) eng_vdwl_thr[tid] += phi;
      if (eflag_atom) eatom_thr[tid][i] += phi;
    }
  }

  // wait until all theads are done with computation
  sync_threads();

  // communicate derivative of embedding function
  // MPI communication only on master thread
#if defined(_OPENMP)
#pragma omp master
#endif
  { comm->forward_comm_pair(this); }

  // wait until master thread is done with communication
  sync_threads();

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

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
	jtype = type[j];
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
	// z2 = phi * r
	// z2p = (phi * r)' = (phi' r) + phi
        // u2 = u
        // u2p = u'
        // w2 = w
        // w2p = w'
	// psip needs both fp[i] and fp[j] terms since r_ij appears in two
	//   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
	//   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

	coeff = rhor_spline[type2rhor[itype][jtype]][m];
	rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
	coeff = rhor_spline[type2rhor[jtype][itype]][m];
	rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
	coeff = z2r_spline[type2z2r[itype][jtype]][m];
	z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
	z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = u2r_spline[type2u2r[itype][jtype]][m];
	u2p = (coeff[0]*p + coeff[1])*p + coeff[2];
	u2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = w2r_spline[type2w2r[itype][jtype]][m];
	w2p = (coeff[0]*p + coeff[1])*p + coeff[2];
	w2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

	recip = 1.0/r;
	phi = z2*recip;
	phip = z2p*recip - phi*recip;
	psip = fp[i]*rhojp + fp[j]*rhoip + phip;
	fpair = -psip*recip;

        delmux = mu[i][0]-mu[j][0];
	delmuy = mu[i][1]-mu[j][1];
	delmuz = mu[i][2]-mu[j][2];
        trdelmu = delmux*delx+delmuy*dely+delmuz*delz;
        sumlamxx = lambda[i][0]+lambda[j][0];
        sumlamyy = lambda[i][1]+lambda[j][1];
        sumlamzz = lambda[i][2]+lambda[j][2];
        sumlamyz = lambda[i][3]+lambda[j][3];
        sumlamxz = lambda[i][4]+lambda[j][4];
        sumlamxy = lambda[i][5]+lambda[j][5];
        tradellam = sumlamxx*delx*delx+sumlamyy*dely*dely+
	  sumlamzz*delz*delz+2.0*sumlamxy*delx*dely+
	  2.0*sumlamxz*delx*delz+2.0*sumlamyz*dely*delz;
        nu = sumlamxx+sumlamyy+sumlamzz;

        adpx = delmux*u2 + trdelmu*u2p*delx*recip + 
	  2.0*w2*(sumlamxx*delx+sumlamxy*dely+sumlamxz*delz) + 
	  w2p*delx*recip*tradellam - 1.0/3.0*nu*(w2p*r+2.0*w2)*delx;
        adpy = delmuy*u2 + trdelmu*u2p*dely*recip + 
	  2.0*w2*(sumlamxy*delx+sumlamyy*dely+sumlamyz*delz) + 
	  w2p*dely*recip*tradellam - 1.0/3.0*nu*(w2p*r+2.0*w2)*dely;
        adpz = delmuz*u2 + trdelmu*u2p*delz*recip + 
	  2.0*w2*(sumlamxz*delx+sumlamyz*dely+sumlamzz*delz) + 
	  w2p*delz*recip*tradellam - 1.0/3.0*nu*(w2p*r+2.0*w2)*delz;
        adpx*=-1.0; adpy*=-1.0; adpz*=-1.0;

        fx = delx*fpair+adpx;
	fy = dely*fpair+adpy;
	fz = delz*fpair+adpz;

	fxtmp += fx;
	fytmp += fy;
	fztmp += fz;
	if (NEWTON_PAIR || j < nlocal) {
	  f[j][0] -= fx;
	  f[j][1] -= fy;
	  f[j][2] -= fz;
	}

	if (EFLAG) evdwl = phi;
	if (EVFLAG) ev_tally_xyz_thr(this,i,j,nlocal,NEWTON_PAIR,evdwl,0.0,
				     fx,fy,fz,delx,dely,delz,tid);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairADPOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairADP::memory_usage();

  return bytes;
}

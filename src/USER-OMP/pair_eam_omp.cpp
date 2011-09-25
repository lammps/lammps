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

#include "pair_eam_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEAMOMP::PairEAMOMP(LAMMPS *lmp) :
  PairEAM(lmp), ThrOMP(lmp, PAIR)
{
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairEAMOMP::compute(int eflag, int vflag)
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
    memory->create(fp,nmax,"pair:fp");
  }

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int ifrom, ito, tid;
    double **f, *rho_t;

    f = loop_setup_thr(atom->f, ifrom, ito, tid, inum, nall, nthreads);
    if (force->newton_pair)
      rho_t = rho + tid*nall;
    else rho_t = rho + tid*atom->nlocal;

    if (evflag) {
      if (eflag) {
	if (force->newton_pair) eval<1,1,1>(f, rho_t, ifrom, ito, tid);
	else eval<1,1,0>(f, rho_t, ifrom, ito, tid);
      } else {
	if (force->newton_pair) eval<1,0,1>(f, rho_t, ifrom, ito, tid);
	else eval<1,0,0>(f, rho_t, ifrom, ito, tid);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(f, rho_t, ifrom, ito, tid);
      else eval<0,0,0>(f, rho_t, ifrom, ito, tid);
    }

    // reduce per thread forces into global force array.
    data_reduce_thr(&(atom->f[0][0]), nall, nthreads, 3, tid);
  } // end of omp parallel region

  // reduce per thread energy and virial, if requested.
  if (evflag) ev_reduce_thr(this);
  if (vflag_fdotr) virial_fdotr_compute();
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairEAMOMP::eval(double **f, double *rho_t,
		      int iifrom, int iito, int tid)
{
  int i,j,ii,jj,m,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
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

  // zero out density 

  if (NEWTON_PAIR) memset(rho_t, 0, nall*sizeof(double));
  else memset(rho_t, 0, nlocal*sizeof(double));

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
	if (NEWTON_PAIR || j < nlocal) {
	  coeff = rhor_spline[type2rhor[itype][jtype]][m];
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
    { comm->reverse_comm_pair(this); }

    // wait until master thread is done with communication
    sync_threads();
  
  } else {
    data_reduce_thr(&(rho[0]), nlocal, comm->nthreads, 1, tid);

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

	recip = 1.0/r;
	phi = z2*recip;
	phip = z2p*recip - phi*recip;
	psip = fp[i]*rhojp + fp[j]*rhoip + phip;
	fpair = -psip*recip;

	fxtmp += delx*fpair;
	fytmp += dely*fpair;
	fztmp += delz*fpair;
	if (NEWTON_PAIR || j < nlocal) {
	  f[j][0] -= delx*fpair;
	  f[j][1] -= dely*fpair;
	  f[j][2] -= delz*fpair;
	}

	if (EFLAG) evdwl = phi;
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

double PairEAMOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairEAM::memory_usage();

  return bytes;
}

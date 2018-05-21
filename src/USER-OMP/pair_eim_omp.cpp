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

#include <cmath>
#include <cstring>

#include "pair_eim_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEIMOMP::PairEIMOMP(LAMMPS *lmp) :
  PairEIM(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairEIMOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
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
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (force->newton_pair)
      thr->init_eim(nall, rho, fp);
    else
      thr->init_eim(atom->nlocal, rho, fp);

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

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairEIMOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,m,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,phip,phi,coul,coulp,recip,psip;
  double *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;


  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  double * const rho_t = thr->get_rho();
  double * const fp_t = thr->get_fp();
  const int tid = thr->get_tid();
  const int nthreads = comm->nthreads;

  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
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
    thr->timer(Timer::PAIR);
    data_reduce_thr(rho, nall, nthreads, 1, tid);

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
    thr->timer(Timer::PAIR);
    data_reduce_thr(rho, nlocal, nthreads, 1, tid);

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
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
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
    thr->timer(Timer::PAIR);
    data_reduce_thr(fp, nall, nthreads, 1, tid);

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
    thr->timer(Timer::PAIR);
    data_reduce_thr(fp, nlocal, nthreads, 1, tid);

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
      e_tally_thr(this, i, i, nlocal, NEWTON_PAIR, phi, 0.0, thr);
    }
  }

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    fxtmp = fytmp = fztmp = 0.0;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
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
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG) evdwl = phi-q0[itype]*q0[jtype]*coul;
        if (EVFLAG) ev_tally_thr(this, i,j,nlocal,NEWTON_PAIR,
                                 evdwl,0.0,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairEIMOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairEIM::memory_usage();

  return bytes;
}

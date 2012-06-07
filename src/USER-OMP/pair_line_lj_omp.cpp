/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "pair_line_lj_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"

#include <string.h>

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLineLJOMP::PairLineLJOMP(LAMMPS *lmp) :
  PairLineLJ(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLineLJOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;
  const int * const line = atom->line;
  const int * const type = atom->type;

  // grow discrete list if necessary and initialize

  if (nall > nmax) {
    nmax = nall;
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->create(dnum,nall,"pair:dnum");
    memory->create(dfirst,nall,"pair:dfirst");
  }
  memset(dnum,0,nall*sizeof(int));
  ndiscrete = 0;

  // need to discretize the system ahead of time
  // until we find a good way to multi-thread it.
  for (int i = 0; i < nall; ++i)
    if (line[i] >= 0)
      if (dnum[i] == 0)
        discretize(i,sigma[type[i]][type[i]]);

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
void PairLineLJOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  int ni,nj,npi,npj,ifirst,jfirst;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,term1,term2,sig,sig3,forcelj;
  double xi[2],xj[2],fi[2],fj[2],dxi,dxj,dyi,dyj,ti,tj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  double * const * const torque = thr->get_torque();
  const int * const line = atom->line;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

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
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      // line/line interactions = NxN particles

      evdwl = 0.0;
      if (line[i] >= 0 && line[j] >= 0) {
        npi = dnum[i];
        ifirst = dfirst[i];
        npj = dnum[j];
        jfirst = dfirst[j];

        fi[0] = fi[1] = fj[0] = fj[1] = ti = tj = 0.0;

        for (ni = 0; ni < npi; ni++) {
          dxi = discrete[ifirst+ni].dx;
          dyi = discrete[ifirst+ni].dy;

          for (nj = 0; nj < npj; nj++) {
            dxj = discrete[jfirst+nj].dx;
            dyj = discrete[jfirst+nj].dy;

            xi[0] = x[i][0] + dxi;
            xi[1] = x[i][1] + dyi;
            xj[0] = x[j][0] + dxj;
            xj[1] = x[j][1] + dyj;

            delx = xi[0] - xj[0];
            dely = xi[1] - xj[1];
            rsq = delx*delx + dely*dely;

            sig = 0.5 * (discrete[ifirst+ni].sigma+discrete[jfirst+nj].sigma);
            sig3 = sig*sig*sig;
            term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
            term1 = 2.0 * term2 * sig3*sig3;
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (term1*r6inv - term2);
            fpair = forcelj*r2inv;

            if (EFLAG) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

            fi[0] += delx*fpair;
            fi[1] += dely*fpair;
            ti += fpair*(dxi*dely - dyi*delx);

            if (NEWTON_PAIR || j < nlocal) {
              fj[0] -= delx*fpair;
              fj[1] -= dely*fpair;
              tj += fpair*(dxj*dely - dyj*delx);
            }
          }
        }

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        torque[i][2] += ti;
        torque[j][2] += tj;

      // line/particle interaction = Nx1 particles
      // convert line into Np particles based on sigma and line length

      } else if (line[i] >= 0) {
        npi = dnum[i];
        ifirst = dfirst[i];

        fi[0] = fi[1] = fj[0] = fj[1] = ti = tj = 0.0;

        for (ni = 0; ni < npi; ni++) {
          dxi = discrete[ifirst+ni].dx;
          dyi = discrete[ifirst+ni].dy;

          xi[0] = x[i][0] + dxi;
          xi[1] = x[i][1] + dyi;
          xj[0] = x[j][0];
          xj[1] = x[j][1];

          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          rsq = delx*delx + dely*dely;

          sig = 0.5 * (discrete[ifirst+ni].sigma+sigma[jtype][jtype]);
          sig3 = sig*sig*sig;
          term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
          term1 = 2.0 * term2 * sig3*sig3;
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (term1*r6inv - term2);
          fpair = forcelj*r2inv;

          if (EFLAG) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

          fi[0] += delx*fpair;
          fi[1] += dely*fpair;
          ti += fpair*(dxi*dely - dyi*delx);

          if (NEWTON_PAIR || j < nlocal) {
            fj[0] -= delx*fpair;
            fj[1] -= dely*fpair;
            tj += fpair*(dxj*dely - dyj*delx);
          }
        }

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        torque[i][2] += ti;
        torque[j][2] += tj;

      // particle/line interaction = Nx1 particles
      // convert line into Np particles based on sigma and line length

      } else if (line[j] >= 0) {
        npj = dnum[j];
        jfirst = dfirst[j];

        fi[0] = fi[1] = fj[0] = fj[1] = ti = tj = 0.0;

        for (nj = 0; nj < npj; nj++) {
          dxj = discrete[jfirst+nj].dx;
          dyj = discrete[jfirst+nj].dy;

          xi[0] = x[i][0];
          xi[1] = x[i][1];
          xj[0] = x[j][0] + dxj;
          xj[1] = x[j][1] + dyj;

          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          rsq = delx*delx + dely*dely;

          sig = 0.5 * (sigma[itype][itype]+discrete[jfirst+nj].sigma);
          sig3 = sig*sig*sig;
          term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
          term1 = 2.0 * term2 * sig3*sig3;
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (term1*r6inv - term2);
          fpair = forcelj*r2inv;

          if (EFLAG) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

          fi[0] += delx*fpair;
          fi[1] += dely*fpair;
          ti += fpair*(dxi*dely - dyi*delx);

          if (NEWTON_PAIR || j < nlocal) {
            fj[0] -= delx*fpair;
            fj[1] -= dely*fpair;
            tj -= fpair*(dxj*dely - dyj*delx);
          }
        }

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        torque[i][2] += ti;
        torque[j][2] += tj;

      // particle/particle interaction = 1x1 particles

      } else {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = forcelj*r2inv;

        if (EFLAG)
          evdwl += r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
      }

      if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                               evdwl,0.0,fpair,delx,dely,delz,thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairLineLJOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLineLJ::memory_usage();

  return bytes;
}

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
#include "pair_tri_lj_omp.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"

#include <string.h>

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairTriLJOMP::PairTriLJOMP(LAMMPS *lmp) :
  PairTriLJ(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairTriLJOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;
  const int * const tri = atom->tri;
  const int * const type = atom->type;
  AtomVecTri::Bonus * const bonus = avec->bonus;

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
  for (int i = 0; i < nall; ++i) {
    double dc1[3],dc2[3],dc3[3],p[3][3];

    if (tri[i] >= 0) {
      if (dnum[i] == 0) {
        MathExtra::quat_to_mat(bonus[tri[i]].quat,p);
        MathExtra::matvec(p,bonus[tri[i]].c1,dc1);
        MathExtra::matvec(p,bonus[tri[i]].c2,dc2);
        MathExtra::matvec(p,bonus[tri[i]].c3,dc3);
        dfirst[i] = ndiscrete;
        discretize(i,sigma[type[i]][type[i]],dc1,dc2,dc3);
        dnum[i] = ndiscrete - dfirst[i];
      }
    }
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
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairTriLJOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  int ni,nj,npi,npj,ifirst,jfirst;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,term1,term2,sig,sig3,forcelj;
  double dxi,dxj,dyi,dyj,dzi,dzj;
  double xi[3],xj[3],fi[3],fj[3],ti[3],tj[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  double * const * const torque = thr->get_torque();
  const int * const tri = atom->tri;
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

      // tri/tri interactions = NxN particles
      // c1,c2,c3 = corner pts of triangle I or J

      evdwl = 0.0;
      if (tri[i] >= 0 && tri[j] >= 0) {
        npi = dnum[i];
        ifirst = dfirst[i];
        npj = dnum[j];
        jfirst = dfirst[j];

        fi[0]=fi[1]=fi[2]=fj[0]=fj[1]=fj[2]=0.0;
        ti[0]=ti[1]=ti[2]=tj[0]=tj[1]=tj[2]=0.0;

        for (ni = 0; ni < npi; ni++) {
          dxi = discrete[ifirst+ni].dx;
          dyi = discrete[ifirst+ni].dy;
          dzi = discrete[ifirst+ni].dz;

          for (nj = 0; nj < npj; nj++) {
            dxj = discrete[jfirst+nj].dx;
            dyj = discrete[jfirst+nj].dy;
            dzj = discrete[jfirst+nj].dz;

            xi[0] = x[i][0] + dxi;
            xi[1] = x[i][1] + dyi;
            xi[2] = x[i][2] + dzi;
            xj[0] = x[j][0] + dxj;
            xj[1] = x[j][1] + dyj;
            xj[2] = x[j][2] + dzj;

            delx = xi[0] - xj[0];
            dely = xi[1] - xj[1];
            delz = xi[2] - xj[2];
            rsq = delx*delx + dely*dely + delz*delz;

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
            fi[2] += delz*fpair;
            ti[0] += fpair*(dyi*delz - dzi*dely);
            ti[1] += fpair*(dzi*delx - dxi*delz);
            ti[2] += fpair*(dxi*dely - dyi*delx);

            if (NEWTON_PAIR || j < nlocal) {
              fj[0] -= delx*fpair;
              fj[1] -= dely*fpair;
              fj[2] -= delz*fpair;
              tj[0] -= fpair*(dyj*delz - dzj*dely);
              tj[1] -= fpair*(dzj*delx - dxj*delz);
              tj[2] -= fpair*(dxj*dely - dyj*delx);
            }
          }
        }

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[i][2] += fi[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        torque[i][0] += ti[0];
        torque[i][1] += ti[1];
        torque[i][2] += ti[2];
        torque[j][0] += tj[0];
        torque[j][1] += tj[1];
        torque[j][2] += tj[2];

        // tri/particle interaction = Nx1 particles
        // c1,c2,c3 = corner pts of triangle I

      } else if (tri[i] >= 0) {
        npi = dnum[i];
        ifirst = dfirst[i];

        fi[0]=fi[1]=fi[2]=fj[0]=fj[1]=fj[2]=0.0;
        ti[0]=ti[1]=ti[2]=tj[0]=tj[1]=tj[2]=0.0;

        for (ni = 0; ni < npi; ni++) {
          dxi = discrete[ifirst+ni].dx;
          dyi = discrete[ifirst+ni].dy;
          dzi = discrete[ifirst+ni].dz;

          xi[0] = x[i][0] + dxi;
          xi[1] = x[i][1] + dyi;
          xi[2] = x[i][2] + dzi;
          xj[0] = x[j][0];
          xj[1] = x[j][1];
          xj[2] = x[j][2];

          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          delz = xi[2] - xj[2];
          rsq = delx*delx + dely*dely + delz*delz;

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
          fi[2] += delz*fpair;
          ti[0] += fpair*(dyi*delz - dzi*dely);
          ti[1] += fpair*(dzi*delx - dxi*delz);
          ti[2] += fpair*(dxi*dely - dyi*delx);

          if (NEWTON_PAIR || j < nlocal) {
            fj[0] -= delx*fpair;
            fj[1] -= dely*fpair;
            fj[2] -= delz*fpair;
            tj[0] -= fpair*(dyj*delz - dzj*dely);
            tj[1] -= fpair*(dzj*delx - dxj*delz);
            tj[2] -= fpair*(dxj*dely - dyj*delx);
          }
        }

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[i][2] += fi[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        torque[i][0] += ti[0];
        torque[i][1] += ti[1];
        torque[i][2] += ti[2];
        torque[j][0] += tj[0];
        torque[j][1] += tj[1];
        torque[j][2] += tj[2];

        // particle/tri interaction = Nx1 particles
        // c1,c2,c3 = corner pts of triangle J

      } else if (tri[j] >= 0) {
        npj = dnum[j];
        jfirst = dfirst[j];

        fi[0]=fi[1]=fi[2]=fj[0]=fj[1]=fj[2]=0.0;
        ti[0]=ti[1]=ti[2]=tj[0]=tj[1]=tj[2]=0.0;

        for (nj = 0; nj < npj; nj++) {
          dxj = discrete[jfirst+nj].dx;
          dyj = discrete[jfirst+nj].dy;
          dzj = discrete[jfirst+nj].dz;

          xi[0] = x[i][0];
          xi[1] = x[i][1];
          xi[2] = x[i][2];
          xj[0] = x[j][0] + dxj;
          xj[1] = x[j][1] + dyj;
          xj[2] = x[j][2] + dzj;

          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          delz = xi[2] - xj[2];
          rsq = delx*delx + dely*dely + delz*delz;

          sig = 0.5 * (sigma[itype][itype]+discrete[jfirst+nj].sigma);
          sig3 = sig*sig*sig;
          term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
          term1 = 2.0 * term2 * sig3*sig3;
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (term1*r6inv - term2);
          fpair = forcelj*r2inv;

          if (EFLAG) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

          if (EFLAG) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

          fi[0] += delx*fpair;
          fi[1] += dely*fpair;
          fi[2] += delz*fpair;
          ti[0] += fpair*(dyi*delz - dzi*dely);
          ti[1] += fpair*(dzi*delx - dxi*delz);
          ti[2] += fpair*(dxi*dely - dyi*delx);

          if (NEWTON_PAIR || j < nlocal) {
            fj[0] -= delx*fpair;
            fj[1] -= dely*fpair;
            fj[2] -= delz*fpair;
            tj[0] -= fpair*(dyj*delz - dzj*dely);
            tj[1] -= fpair*(dzj*delx - dxj*delz);
            tj[2] -= fpair*(dxj*dely - dyj*delx);
          }
        }

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[i][2] += fi[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        torque[i][0] += ti[0];
        torque[i][1] += ti[1];
        torque[i][2] += ti[2];
        torque[j][0] += tj[0];
        torque[j][1] += tj[1];
        torque[j][2] += tj[2];

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

double PairTriLJOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairTriLJ::memory_usage();

  return bytes;
}

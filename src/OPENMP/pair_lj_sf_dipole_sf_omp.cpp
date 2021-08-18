// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_lj_sf_dipole_sf_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJSFDipoleSFOMP::PairLJSFDipoleSFOMP(LAMMPS *lmp) :
  PairLJSFDipoleSF(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJSFDipoleSFOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

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
void PairLJSFDipoleSFOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,qtmp,delx,dely,delz,evdwl,ecoul;
  double rsq,rinv,r2inv,r6inv,r3inv,r5inv,fx,fy,fz;
  double forcecoulx,forcecouly,forcecoulz,crossx,crossy,crossz;
  double tixcoul,tiycoul,tizcoul,tjxcoul,tjycoul,tjzcoul;
  double fq,pdotp,pidotr,pjdotr,pre1,pre2,pre3,pre4;
  double forcelj,factor_coul,factor_lj;
  double presf,afac,bfac,pqfac,qpfac,forceljcut,forceljsf;
  double aforcecoulx,aforcecouly,aforcecoulz;
  double bforcecoulx,bforcecouly,bforcecoulz;
  double rcutlj2inv, rcutcoul2inv,rcutlj6inv;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  double * const * const torque = thr->get_torque();
  const double * _noalias const q = atom->q;
  const dbl4_t * _noalias const mu = (dbl4_t *) atom->mu[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const double * _noalias const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  double fxtmp,fytmp,fztmp,t1tmp,t2tmp,t3tmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    qtmp = q[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=t1tmp=t2tmp=t3tmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        rinv = sqrt(r2inv);

        // atom can have both a charge and dipole
        // i,j = charge-charge, dipole-dipole, dipole-charge, or charge-dipole

        forcecoulx = forcecouly = forcecoulz = 0.0;
        tixcoul = tiycoul = tizcoul = 0.0;
        tjxcoul = tjycoul = tjzcoul = 0.0;

        if (rsq < cut_coulsq[itype][jtype]) {

          rcutcoul2inv=1.0/cut_coulsq[itype][jtype];
          if (qtmp != 0.0 && q[j] != 0.0) {
            pre1 = qtmp*q[j]*rinv*(r2inv-rcutcoul2inv);

            forcecoulx += pre1*delx;
            forcecouly += pre1*dely;
            forcecoulz += pre1*delz;
          }

          if (mu[i].w > 0.0 && mu[j].w > 0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;

            pdotp = mu[i].x*mu[j].x + mu[i].y*mu[j].y + mu[i].z*mu[j].z;
            pidotr = mu[i].x*delx + mu[i].y*dely + mu[i].z*delz;
            pjdotr = mu[j].x*delx + mu[j].y*dely + mu[j].z*delz;

            afac = 1.0 - rsq*rsq * rcutcoul2inv*rcutcoul2inv;
            pre1 = afac * ( pdotp - 3.0 * r2inv * pidotr * pjdotr );
            aforcecoulx = pre1*delx;
            aforcecouly = pre1*dely;
            aforcecoulz = pre1*delz;

            bfac = 1.0 - 4.0*rsq*sqrt(rsq*rcutcoul2inv)*rcutcoul2inv +
              3.0*rsq*rsq*rcutcoul2inv*rcutcoul2inv;
            presf = 2.0 * r2inv * pidotr * pjdotr;
            bforcecoulx = bfac * (pjdotr*mu[i].x+pidotr*mu[j].x-presf*delx);
            bforcecouly = bfac * (pjdotr*mu[i].y+pidotr*mu[j].y-presf*dely);
            bforcecoulz = bfac * (pjdotr*mu[i].z+pidotr*mu[j].z-presf*delz);

            forcecoulx += 3.0 * r5inv * ( aforcecoulx + bforcecoulx );
            forcecouly += 3.0 * r5inv * ( aforcecouly + bforcecouly );
            forcecoulz += 3.0 * r5inv * ( aforcecoulz + bforcecoulz );

            pre2 = 3.0 * bfac * r5inv * pjdotr;
            pre3 = 3.0 * bfac * r5inv * pidotr;
            pre4 = -bfac * r3inv;

            crossx = pre4 * (mu[i].y*mu[j].z - mu[i].z*mu[j].y);
            crossy = pre4 * (mu[i].z*mu[j].x - mu[i].x*mu[j].z);
            crossz = pre4 * (mu[i].x*mu[j].y - mu[i].y*mu[j].x);

            tixcoul += crossx + pre2 * (mu[i].y*delz - mu[i].z*dely);
            tiycoul += crossy + pre2 * (mu[i].z*delx - mu[i].x*delz);
            tizcoul += crossz + pre2 * (mu[i].x*dely - mu[i].y*delx);
            tjxcoul += -crossx + pre3 * (mu[j].y*delz - mu[j].z*dely);
            tjycoul += -crossy + pre3 * (mu[j].z*delx - mu[j].x*delz);
            tjzcoul += -crossz + pre3 * (mu[j].x*dely - mu[j].y*delx);
          }

          if (mu[i].w > 0.0 && q[j] != 0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pidotr = mu[i].x*delx + mu[i].y*dely + mu[i].z*delz;
            pre1 = 3.0 * q[j] * r5inv * pidotr * (1-rsq*rcutcoul2inv);
            pqfac = 1.0 - 3.0*rsq*rcutcoul2inv +
              2.0*rsq*sqrt(rsq*rcutcoul2inv)*rcutcoul2inv;
            pre2 = q[j] * r3inv * pqfac;

            forcecoulx += pre2*mu[i].x - pre1*delx;
            forcecouly += pre2*mu[i].y - pre1*dely;
            forcecoulz += pre2*mu[i].z - pre1*delz;
            tixcoul += pre2 * (mu[i].y*delz - mu[i].z*dely);
            tiycoul += pre2 * (mu[i].z*delx - mu[i].x*delz);
            tizcoul += pre2 * (mu[i].x*dely - mu[i].y*delx);
          }

          if (mu[j].w > 0.0 && qtmp != 0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pjdotr = mu[j].x*delx + mu[j].y*dely + mu[j].z*delz;
            pre1 = 3.0 * qtmp * r5inv * pjdotr * (1-rsq*rcutcoul2inv);
            qpfac = 1.0 - 3.0*rsq*rcutcoul2inv +
              2.0*rsq*sqrt(rsq*rcutcoul2inv)*rcutcoul2inv;
            pre2 = qtmp * r3inv * qpfac;

            forcecoulx += pre1*delx - pre2*mu[j].x;
            forcecouly += pre1*dely - pre2*mu[j].y;
            forcecoulz += pre1*delz - pre2*mu[j].z;
            tjxcoul += -pre2 * (mu[j].y*delz - mu[j].z*dely);
            tjycoul += -pre2 * (mu[j].z*delx - mu[j].x*delz);
            tjzcoul += -pre2 * (mu[j].x*dely - mu[j].y*delx);
          }
        }

        // LJ interaction

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          forceljcut = r6inv*(lj1[itype][jtype]*r6inv-lj2[itype][jtype])*r2inv;

          rcutlj2inv = 1.0 / cut_ljsq[itype][jtype];
          rcutlj6inv = rcutlj2inv * rcutlj2inv * rcutlj2inv;
          forceljsf = (lj1[itype][jtype]*rcutlj6inv - lj2[itype][jtype]) *
            rcutlj6inv * rcutlj2inv;

          forcelj = factor_lj * (forceljcut - forceljsf);
        } else forcelj = 0.0;

        // total force

        fq = factor_coul*qqrd2e*scale[itype][jtype];
        fx = fq*forcecoulx + delx*forcelj;
        fy = fq*forcecouly + dely*forcelj;
        fz = fq*forcecoulz + delz*forcelj;

        // force & torque accumulation

        fxtmp += fx;
        fytmp += fy;
        fztmp += fz;
        t1tmp += fq*tixcoul;
        t2tmp += fq*tiycoul;
        t3tmp += fq*tizcoul;

        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= fx;
          f[j].y -= fy;
          f[j].z -= fz;
          torque[j][0] += fq*tjxcoul;
          torque[j][1] += fq*tjycoul;
          torque[j][2] += fq*tjzcoul;
        }

        if (EFLAG) {
          if (rsq < cut_coulsq[itype][jtype]) {
            ecoul = (1.0-sqrt(rsq/cut_coulsq[itype][jtype]));
            ecoul *= ecoul;
            ecoul *= qtmp * q[j] * rinv;
            if (mu[i].w > 0.0 && mu[j].w > 0.0)
              ecoul += bfac * (r3inv*pdotp - 3.0*r5inv*pidotr*pjdotr);
            if (mu[i].w > 0.0 && q[j] != 0.0)
              ecoul += -q[j] * r3inv * pqfac * pidotr;
            if (mu[j].w > 0.0 && qtmp != 0.0)
              ecoul += qtmp * r3inv * qpfac * pjdotr;
            ecoul *= factor_coul*qqrd2e*scale[itype][jtype];
          } else ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype])+
              rcutlj6inv*(6*lj3[itype][jtype]*rcutlj6inv-3*lj4[itype][jtype])*
              rsq*rcutlj2inv+
              rcutlj6inv*(-7*lj3[itype][jtype]*rcutlj6inv+4*lj4[itype][jtype]);
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (EVFLAG) ev_tally_xyz_thr(this,i,j,nlocal,NEWTON_PAIR,
                                     evdwl,ecoul,fx,fy,fz,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
    torque[i][0] += t1tmp;
    torque[i][1] += t2tmp;
    torque[i][2] += t3tmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJSFDipoleSFOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJSFDipoleSF::memory_usage();

  return bytes;
}

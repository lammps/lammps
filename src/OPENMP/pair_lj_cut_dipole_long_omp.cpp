// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_lj_cut_dipole_long_omp.h"

#include "atom.h"
#include "comm.h"
#include "ewald_const.h"
#include "force.h"
#include "math_const.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

PairLJCutDipoleLongOMP::PairLJCutDipoleLongOMP(LAMMPS *lmp) :
  PairLJCutDipoleLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJCutDipoleLongOMP::compute(int eflag, int vflag)
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

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJCutDipoleLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r,rinv,r2inv,r6inv;
  double forcecoulx,forcecouly,forcecoulz,fforce;
  double tixcoul,tiycoul,tizcoul,tjxcoul,tjycoul,tjzcoul;
  double fx,fy,fz,fdx,fdy,fdz,fax,fay,faz;
  double pdotp,pidotr,pjdotr,pre1,pre2,pre3;
  double grij,expm2,t,erfc;
  double g0,g1,g2,b0,b1,b2,b3,d0,d1,d2,d3;
  double zdix,zdiy,zdiz,zdjx,zdjy,zdjz,zaix,zaiy,zaiz,zajx,zajy,zajz;
  double g0b1_g1b2_g2b3,g0d1_g1d2_g2d3;
  double forcelj,factor_coul,factor_lj,facm1;
  double evdwl,ecoul;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const double * const * const mu = atom->mu;
  double * const * const torque = thr->get_torque();
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  pre1 = 2.0 * g_ewald / MY_PIS;
  pre2 = 4.0 * pow(g_ewald,3.0) / MY_PIS;
  pre3 = 8.0 * pow(g_ewald,5.0) / MY_PIS;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        rinv = sqrt(r2inv);

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

          pdotp = mu[i][0]*mu[j][0] + mu[i][1]*mu[j][1] + mu[i][2]*mu[j][2];
          pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
          pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;

          g0 = qtmp*q[j];
          g1 = qtmp*pjdotr - q[j]*pidotr + pdotp;
          g2 = -pidotr*pjdotr;

          if (factor_coul > 0.0) {
            b0 = erfc * rinv;
            b1 = (b0 + pre1*expm2) * r2inv;
            b2 = (3.0*b1 + pre2*expm2) * r2inv;
            b3 = (5.0*b2 + pre3*expm2) * r2inv;

            g0b1_g1b2_g2b3 = g0*b1 + g1*b2 + g2*b3;
            fdx = delx * g0b1_g1b2_g2b3 -
              b1 * (qtmp*mu[j][0] - q[j]*mu[i][0]) +
              b2 * (pjdotr*mu[i][0] + pidotr*mu[j][0]);
            fdy = dely * g0b1_g1b2_g2b3 -
              b1 * (qtmp*mu[j][1] - q[j]*mu[i][1]) +
              b2 * (pjdotr*mu[i][1] + pidotr*mu[j][1]);
            fdz = delz * g0b1_g1b2_g2b3 -
              b1 * (qtmp*mu[j][2] - q[j]*mu[i][2]) +
              b2 * (pjdotr*mu[i][2] + pidotr*mu[j][2]);

            zdix = delx * (q[j]*b1 + b2*pjdotr) - b1*mu[j][0];
            zdiy = dely * (q[j]*b1 + b2*pjdotr) - b1*mu[j][1];
            zdiz = delz * (q[j]*b1 + b2*pjdotr) - b1*mu[j][2];
            zdjx = delx * (-qtmp*b1 + b2*pidotr) - b1*mu[i][0];
            zdjy = dely * (-qtmp*b1 + b2*pidotr) - b1*mu[i][1];
            zdjz = delz * (-qtmp*b1 + b2*pidotr) - b1*mu[i][2];

            if (factor_coul < 1.0) {
              fdx *= factor_coul;
              fdy *= factor_coul;
              fdz *= factor_coul;
              zdix *= factor_coul;
              zdiy *= factor_coul;
              zdiz *= factor_coul;
              zdjx *= factor_coul;
              zdjy *= factor_coul;
              zdjz *= factor_coul;
            }
          } else {
            fdx = fdy = fdz = 0.0;
            zdix = zdiy = zdiz = 0.0;
            zdjx = zdjy = zdjz = 0.0;
          }

          if (factor_coul < 1.0) {
            d0 = (erfc - 1.0) * rinv;
            d1 = (d0 + pre1*expm2) * r2inv;
            d2 = (3.0*d1 + pre2*expm2) * r2inv;
            d3 = (5.0*d2 + pre3*expm2) * r2inv;

            g0d1_g1d2_g2d3 = g0*d1 + g1*d2 + g2*d3;
            fax = delx * g0d1_g1d2_g2d3 -
              d1 * (qtmp*mu[j][0] - q[j]*mu[i][0]) +
              d2 * (pjdotr*mu[i][0] + pidotr*mu[j][0]);
            fay = dely * g0d1_g1d2_g2d3 -
              d1 * (qtmp*mu[j][1] - q[j]*mu[i][1]) +
              d2 * (pjdotr*mu[i][1] + pidotr*mu[j][1]);
            faz = delz * g0d1_g1d2_g2d3 -
              d1 * (qtmp*mu[j][2] - q[j]*mu[i][2]) +
              d2 * (pjdotr*mu[i][2] + pidotr*mu[j][2]);

            zaix = delx * (q[j]*d1 + d2*pjdotr) - d1*mu[j][0];
            zaiy = dely * (q[j]*d1 + d2*pjdotr) - d1*mu[j][1];
            zaiz = delz * (q[j]*d1 + d2*pjdotr) - d1*mu[j][2];
            zajx = delx * (-qtmp*d1 + d2*pidotr) - d1*mu[i][0];
            zajy = dely * (-qtmp*d1 + d2*pidotr) - d1*mu[i][1];
            zajz = delz * (-qtmp*d1 + d2*pidotr) - d1*mu[i][2];

            if (factor_coul > 0.0) {
              facm1 = 1.0 - factor_coul;
              fax *= facm1;
              fay *= facm1;
              faz *= facm1;
              zaix *= facm1;
              zaiy *= facm1;
              zaiz *= facm1;
              zajx *= facm1;
              zajy *= facm1;
              zajz *= facm1;
            }
          } else {
            fax = fay = faz = 0.0;
            zaix = zaiy = zaiz = 0.0;
            zajx = zajy = zajz = 0.0;
          }

          forcecoulx = fdx + fax;
          forcecouly = fdy + fay;
          forcecoulz = fdz + faz;

          tixcoul = mu[i][1]*(zdiz + zaiz) - mu[i][2]*(zdiy + zaiy);
          tiycoul = mu[i][2]*(zdix + zaix) - mu[i][0]*(zdiz + zaiz);
          tizcoul = mu[i][0]*(zdiy + zaiy) - mu[i][1]*(zdix + zaix);
          tjxcoul = mu[j][1]*(zdjz + zajz) - mu[j][2]*(zdjy + zajy);
          tjycoul = mu[j][2]*(zdjx + zajx) - mu[j][0]*(zdjz + zajz);
          tjzcoul = mu[j][0]*(zdjy + zajy) - mu[j][1]*(zdjx + zajx);

        } else {
          forcecoulx = forcecouly = forcecoulz = 0.0;
          tixcoul = tiycoul = tizcoul = 0.0;
          tjxcoul = tjycoul = tjzcoul = 0.0;
        }

        // LJ interaction

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          fforce = factor_lj * forcelj*r2inv;
        } else fforce = 0.0;

        // total force

        fx = qqrd2e*forcecoulx + delx*fforce;
        fy = qqrd2e*forcecouly + dely*fforce;
        fz = qqrd2e*forcecoulz + delz*fforce;

        // force & torque accumulation

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        torque[i][0] += qqrd2e*tixcoul;
        torque[i][1] += qqrd2e*tiycoul;
        torque[i][2] += qqrd2e*tizcoul;

        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] += qqrd2e*tjxcoul;
          torque[j][1] += qqrd2e*tjycoul;
          torque[j][2] += qqrd2e*tjzcoul;
        }

        if (EFLAG) {
          if (rsq < cut_coulsq && factor_coul > 0.0) {
            ecoul = qqrd2e*(b0*g0 + b1*g1 + b2*g2);
            if (factor_coul < 1.0) {
              ecoul *= factor_coul;
              ecoul += (1-factor_coul) * qqrd2e * (d0*g0 + d1*g1 + d2*g2);
            }
          } else ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
              offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (EVFLAG) ev_tally_xyz_thr(this,i,j,nlocal,NEWTON_PAIR,
                                     evdwl,ecoul,fx,fy,fz,delx,dely,delz,thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCutDipoleLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJCutDipoleLong::memory_usage();

  return bytes;
}

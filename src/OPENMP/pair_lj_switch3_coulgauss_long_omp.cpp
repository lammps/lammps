// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Steven Vandenbrande
                        (OMP threading by Axel Kohlmeyer, Temple U)
------------------------------------------------------------------------- */

#include "omp_compat.h"
#include "pair_lj_switch3_coulgauss_long_omp.h"

#include "atom.h"
#include "comm.h"
#include "ewald_const.h"
#include "force.h"
#include "math_const.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

PairLJSwitch3CoulGaussLongOMP::PairLJSwitch3CoulGaussLongOMP(LAMMPS *lmp) :
  PairLJSwitch3CoulGaussLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJSwitch3CoulGaussLongOMP::compute(int eflag, int vflag)
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
void PairLJSwitch3CoulGaussLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double r,rsq,r2inv,r6inv,forcecoul,forcecoul2,forcelj,factor_coul,factor_lj,tr,ftr,trx;
  double grij,expm2,prefactor,prefactor2,t,erfc1,erfc2,rrij;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const double * _noalias const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

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
        forcecoul = forcecoul2 = forcelj = 0.0;
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            erfc1 = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc1 + EWALD_F*grij*expm2);
            if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        }

        if (rsq < cut_ljsq[itype][jtype]) {
          // Lennard-Jones potential
          r = sqrt(rsq);
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv*(12.0*lj3[itype][jtype]*r6inv - 6.0*lj4[itype][jtype]);
          // Correction for Gaussian radii
          if (lj2[itype][jtype]==0.0) {
            erfc2 = 0.0;
            prefactor2 = 0.0;
          } else {
            rrij = lj2[itype][jtype]*r;
            double expn2 = exp(-rrij*rrij);
            erfc2 = erfc(rrij);
            prefactor2 = -qqrd2e*qtmp*q[j]/r;
            forcecoul2 = prefactor2*(erfc2+EWALD_F*rrij*expn2);
          }
          // VDW energy must be computed unconditionally: it is needed for the
          // switching/truncation force term even when energy is not requested
          evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
              offset[itype][jtype];
        } else evdwl = 0.0;

        if (EFLAG) {
          if (rsq < cut_coulsq) {
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor*erfc1;
            else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype])
            ecoul += prefactor2*erfc2*factor_coul;
        }

        // Truncation, see Yaff Switch3
        if (truncw>0) {
          if (rsq < cut_ljsq[itype][jtype]) {
            if (r>cut_lj[itype][jtype]-truncw) {
              trx = (cut_lj[itype][jtype]-r)*truncwi;
              tr = trx*trx*(3.0-2.0*trx);
              ftr = 6.0*trx*(1.0-trx)*r*truncwi;
              forcelj = forcelj*tr + evdwl*ftr;
              evdwl *= tr;
            }
          }
        }

        fpair = (forcecoul + factor_coul*forcecoul2 + factor_lj*forcelj) * r2inv;
        evdwl *= factor_lj;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EVFLAG) ev_tally_thr(this, i,j,nlocal,NEWTON_PAIR,
                                 evdwl,ecoul,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJSwitch3CoulGaussLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJSwitch3CoulGaussLong::memory_usage();
  return bytes;
}

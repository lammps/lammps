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
   Contributing authors: Hendrik Heenen (Technical University of Munich)
                         and Rochus Schmid (Ruhr-Universitaet Bochum)
                         (OMP threading by Axel Kohlmeyer, Temple U)
------------------------------------------------------------------------- */

#include "omp_compat.h"
#include "pair_buck6d_coul_gauss_dsf_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairBuck6dCoulGaussDSFOMP::PairBuck6dCoulGaussDSFOMP(LAMMPS *lmp) :
  PairBuck6dCoulGaussDSF(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSFOMP::compute(int eflag, int vflag)
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
void PairBuck6dCoulGaussDSFOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r,rsq,r2inv,r6inv,r14inv,rexp,forcecoul,forcebuck6d,factor_coul,factor_lj;
  double term1,term2,term3,term4,term5;
  double rcu,rqu,sme,smf,ebuck6d;
  double prefactor,erfcc,erfcd,arg;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  const double * _noalias const special_coul = force->special_coul;
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
        r2inv = 1.0/rsq;
        r = sqrt(rsq);

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          r14inv = r6inv*r6inv*r2inv;
          rexp = exp(-r*buck6d2[itype][jtype]);
          term1 = buck6d3[itype][jtype]*r6inv;
          term2 = buck6d4[itype][jtype]*r14inv;
          term3 = term2*term2;
          term4 = 1.0/(1.0 + term2);
          term5 = 1.0/(1.0 + 2.0*term2 + term3);
          forcebuck6d = buck6d1[itype][jtype]*buck6d2[itype][jtype]*r*rexp;
          forcebuck6d -= term1*(6.0*term4 - term5*14.0*term2);
          ebuck6d = buck6d1[itype][jtype]*rexp - term1*term4;

          // smoothing term
          if (rsq > rsmooth_sq[itype][jtype]) {
            rcu = r*rsq;
            rqu = rsq*rsq;
            sme = c5[itype][jtype]*rqu*r + c4[itype][jtype]*rqu + c3[itype][jtype]*rcu +
                  c2[itype][jtype]*rsq + c1[itype][jtype]*r + c0[itype][jtype];
            smf = 5.0*c5[itype][jtype]*rqu + 4.0*c4[itype][jtype]*rcu +
                  3.0*c3[itype][jtype]*rsq + 2.0*c2[itype][jtype]*r + c1[itype][jtype];
            forcebuck6d = forcebuck6d*sme - ebuck6d*smf*r;
            ebuck6d *= sme;
          }
        } else forcebuck6d = 0.0;

        if (rsq < cut_coulsq) {
          prefactor = qqrd2e*qtmp*q[j]/r;

          arg = alpha_ij[itype][jtype]*r;
          erfcd = MathSpecial::expmsq(arg);
          erfcc = 1 - (MathSpecial::my_erfcx(arg) * erfcd);

          forcecoul = prefactor * ((erfcc/r) - (2.0/MY_PIS*alpha_ij[itype][jtype]*erfcd) +
                                                r*f_shift_ij[itype][jtype]) * r;

          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else forcecoul = 0.0;

        fpair = (forcecoul + factor_lj*forcebuck6d) * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG) {
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = ebuck6d - offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;

          if (rsq < cut_coulsq) {
            ecoul = prefactor * (erfcc - r*e_shift_ij[itype][jtype] -
                                 rsq*f_shift_ij[itype][jtype]);
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
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

double PairBuck6dCoulGaussDSFOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBuck6dCoulGaussDSF::memory_usage();
  return bytes;
}

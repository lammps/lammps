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
   This style is a simplified re-implementation of the CG/CMM pair style
------------------------------------------------------------------------- */

#include "pair_lj_sdk_coul_msm.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "kspace.h"
#include "neigh_list.h"
#include "error.h"

#include "lj_sdk_common.h"

using namespace LAMMPS_NS;
using namespace LJSDKParms;

/* ---------------------------------------------------------------------- */

PairLJSDKCoulMSM::PairLJSDKCoulMSM(LAMMPS *lmp) : PairLJSDKCoulLong(lmp)
{
  ewaldflag = pppmflag = 0;
  msmflag = 1;
  respa_enable = 0;
  ftable = NULL;
}

/* ---------------------------------------------------------------------- */

void PairLJSDKCoulMSM::compute(int eflag, int vflag)
{
  if (force->kspace->scalar_pressure_flag)
    error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' with Pair style");

  ev_init(eflag,vflag);

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) eval_msm<1,1,1>();
      else eval_msm<1,1,0>();
    } else {
      if (force->newton_pair) eval_msm<1,0,1>();
      else eval_msm<1,0,0>();
    }
  } else {
    if (force->newton_pair) eval_msm<0,0,1>();
    else eval_msm<0,0,0>();
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJSDKCoulMSM::eval_msm()
{
  int i,ii,j,jj,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double r,rsq,r2inv,forcecoul,forcelj,factor_coul,factor_lj;
  double egamma,fgamma,prefactor;

  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  double fxtmp,fytmp,fztmp;

  const int inum = list->inum;
  const int * const ilist = list->ilist;
  const int * const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp=fytmp=fztmp=0.0;

    const int itype = type[i];
    const int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      forcecoul = forcelj = evdwl = ecoul = 0.0;

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
        const int ljt = lj_type[itype][jtype];

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            prefactor = qqrd2e * qtmp*q[j]/r;
            fgamma = 1.0 + (rsq/cut_coulsq)*force->kspace->dgamma(r/cut_coul);
            forcecoul = prefactor * fgamma;
            if (EFLAG) {
              egamma = 1.0 - (r/cut_coul)*force->kspace->gamma(r/cut_coul);
              ecoul = prefactor*egamma;
            }
            if (factor_coul < 1.0) {
              forcecoul -= (1.0-factor_coul)*prefactor;
              if (EFLAG) ecoul -= (1.0-factor_coul)*prefactor;
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (EFLAG) ecoul = qtmp*q[j] *
              (etable[itable] + fraction*detable[itable]);
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
              if (EFLAG) ecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {

          if (ljt == LJ12_4) {
            const double r4inv=r2inv*r2inv;
            forcelj = r4inv*(lj1[itype][jtype]*r4inv*r4inv
                             - lj2[itype][jtype]);

            if (EFLAG)
              evdwl = r4inv*(lj3[itype][jtype]*r4inv*r4inv
                             - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ9_6) {
            const double r3inv = r2inv*sqrt(r2inv);
            const double r6inv = r3inv*r3inv;
            forcelj = r6inv*(lj1[itype][jtype]*r3inv
                             - lj2[itype][jtype]);
            if (EFLAG)
              evdwl = r6inv*(lj3[itype][jtype]*r3inv
                             - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ12_6) {
            const double r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv*(lj1[itype][jtype]*r6inv
                             - lj2[itype][jtype]);
            if (EFLAG)
              evdwl = r6inv*(lj3[itype][jtype]*r6inv
                             - lj4[itype][jtype]) - offset[itype][jtype];
          }
          forcelj *= factor_lj;
          if (EFLAG) evdwl *= factor_lj;
        }

        fpair = (forcecoul + forcelj) * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJSDKCoulMSM::single(int i, int j, int itype, int jtype,
                                 double rsq,
                                 double factor_coul, double factor_lj,
                                 double &fforce)
{
  double r2inv,r,fgamma,egamma,prefactor;
  double fraction,table,forcecoul,forcelj,phicoul,philj;
  int itable;

  forcecoul = forcelj = phicoul = philj = 0.0;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq) {
    if (!ncoultablebits || rsq <= tabinnersq) {
      r = sqrt(rsq);
      prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
      egamma = 1.0 - (r/cut_coul)*force->kspace->gamma(r/cut_coul);
      fgamma = 1.0 + (rsq/cut_coulsq)*force->kspace->dgamma(r/cut_coul);
      forcecoul = prefactor * fgamma;
      phicoul = prefactor * egamma;
      if (factor_coul < 1.0) {
        forcecoul -= (1.0-factor_coul)*prefactor;
        phicoul -= (1.0-factor_coul)*prefactor;
      }
    } else {
      union_int_float_t rsq_lookup_single;
      rsq_lookup_single.f = rsq;
      itable = rsq_lookup_single.i & ncoulmask;
      itable >>= ncoulshiftbits;
      fraction = (rsq_lookup_single.f - rtable[itable]) * drtable[itable];
      table = ftable[itable] + fraction*dftable[itable];
      forcecoul = atom->q[i]*atom->q[j] * table;
      table = etable[itable] + fraction*detable[itable];
      phicoul = atom->q[i]*atom->q[j] * table;
      if (factor_coul < 1.0) {
        table = ctable[itable] + fraction*dctable[itable];
        prefactor = atom->q[i]*atom->q[j] * table;
        forcecoul -= (1.0-factor_coul)*prefactor;
        phicoul -= (1.0-factor_coul)*prefactor;
      }
    }
  }

  if (rsq < cut_ljsq[itype][jtype]) {
    const int ljt = lj_type[itype][jtype];
    const double ljpow1 = lj_pow1[ljt];
    const double ljpow2 = lj_pow2[ljt];
    const double ljpref = lj_prefact[ljt];

    const double ratio = sigma[itype][jtype]/sqrt(rsq);
    const double eps = epsilon[itype][jtype];

    forcelj = factor_lj * ljpref*eps * (ljpow1*pow(ratio,ljpow1)
                          - ljpow2*pow(ratio,ljpow2))/rsq;
    philj = factor_lj * (ljpref*eps * (pow(ratio,ljpow1) - pow(ratio,ljpow2))
                         - offset[itype][jtype]);
  }

  fforce = (forcecoul + forcelj) * r2inv;

  return phicoul + philj;
}

/* ---------------------------------------------------------------------- */

void *PairLJSDKCoulMSM::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"lj_type") == 0) return (void *) lj_type;
  if (strcmp(str,"lj1") == 0) return (void *) lj1;
  if (strcmp(str,"lj2") == 0) return (void *) lj2;
  if (strcmp(str,"lj3") == 0) return (void *) lj3;
  if (strcmp(str,"lj4") == 0) return (void *) lj4;
  if (strcmp(str,"rminsq") == 0) return (void *) rminsq;
  if (strcmp(str,"emin") == 0) return (void *) emin;

  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  if (strcmp(str,"cut_msm") == 0) return (void *) &cut_coul;
  return NULL;
}


// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Stan Moore (SNL), Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "pair_coul_msm.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "kspace.h"
#include "neigh_list.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulMSM::PairCoulMSM(LAMMPS *lmp) : PairCoulLong(lmp)
{
  ewaldflag = pppmflag = 0;
  msmflag = 1;
}

/* ---------------------------------------------------------------------- */

void PairCoulMSM::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double fraction,table;
  double r,r2inv,forcecoul,factor_coul;
  double egamma,fgamma,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  if (force->kspace->scalar_pressure_flag && vflag) {
    if (vflag > 2)
      error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' "
        "to obtain per-atom virial with kspace_style MSM");

    // must switch on global energy computation if not already on

    if (eflag == 0 || eflag == 2) {
      eflag++;
    }
  }

  ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
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
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {
        r2inv = 1.0/rsq;
        if (!ncoultablebits || rsq <= tabinnersq) {
          r = sqrt(rsq);
          prefactor = qqrd2e * scale[itype][jtype] * qtmp*q[j]/r;
          egamma = 1.0 - (r/cut_coul)*force->kspace->gamma(r/cut_coul);
          fgamma = 1.0 + (rsq/cut_coulsq)*force->kspace->dgamma(r/cut_coul);
          forcecoul = prefactor * fgamma;
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else {
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & ncoulmask;
          itable >>= ncoulshiftbits;
          fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
          table = ftable[itable] + fraction*dftable[itable];
          forcecoul = scale[itype][jtype] * qtmp*q[j] * table;
          if (factor_coul < 1.0) {
            table = ctable[itable] + fraction*dctable[itable];
            prefactor = scale[itype][jtype] * qtmp*q[j] * table;
            forcecoul -= (1.0-factor_coul)*prefactor;
          }
        }

        fpair = forcecoul * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (!ncoultablebits || rsq <= tabinnersq)
            ecoul = prefactor*egamma;
          else {
            table = etable[itable] + fraction*detable[itable];
            ecoul = scale[itype][jtype] * qtmp*q[j] * table;
          }
          if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
        }

        if (force->kspace->scalar_pressure_flag)
          fpair = 0.0;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             0.0,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr && !force->kspace->scalar_pressure_flag)
    virial_fdotr_compute();
  if (force->kspace->scalar_pressure_flag && vflag)
    for (i = 0; i < 3; i++) virial[i] += force->pair->eng_coul/3.0;
}

/* ---------------------------------------------------------------------- */

double PairCoulMSM::single(int i, int j, int /*itype*/, int /*jtype*/,
                            double rsq,
                            double factor_coul, double /*factor_lj*/,
                            double &fforce)
{
  double r2inv,r,egamma,fgamma,prefactor;
  double fraction,table,forcecoul,phicoul;
  int itable;

  r2inv = 1.0/rsq;
  if (!ncoultablebits || rsq <= tabinnersq) {
    r = sqrt(rsq);
    egamma = 1.0 - (r/cut_coul)*force->kspace->gamma(r/cut_coul);
    fgamma = 1.0 + (rsq/cut_coulsq)*force->kspace->dgamma(r/cut_coul);
    prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
    forcecoul = prefactor * fgamma;

    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
  } else {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    itable = rsq_lookup.i & ncoulmask;
    itable >>= ncoulshiftbits;
    fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
    table = ftable[itable] + fraction*dftable[itable];
    forcecoul = atom->q[i]*atom->q[j] * table;
    if (factor_coul < 1.0) {
      table = ctable[itable] + fraction*dctable[itable];
      prefactor = atom->q[i]*atom->q[j] * table;
      forcecoul -= (1.0-factor_coul)*prefactor;
    }
  }
  fforce = forcecoul * r2inv;

  if (!ncoultablebits || rsq <= tabinnersq)
    phicoul = prefactor*egamma;
  else {
    table = etable[itable] + fraction*detable[itable];
    phicoul = atom->q[i]*atom->q[j] * table;
  }
  if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;

  return phicoul;
}

/* ---------------------------------------------------------------------- */

void *PairCoulMSM::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  if (strcmp(str,"scale") == 0) {
    dim = 2;
    return (void *) scale;
  }
  return nullptr;
}

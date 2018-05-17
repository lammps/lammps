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

#include "bond_quartic_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "domain.h"
#include "pair.h"

#include <cmath>

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondQuarticOMP::BondQuarticOMP(class LAMMPS *lmp)
  : BondQuartic(lmp), ThrOMP(lmp,THR_BOND)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void BondQuarticOMP::compute(int eflag, int vflag)
{

  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = 0;

  // insure pair->ev_tally() will use 1-4 virial contribution

  if (vflag_global == 2)
    force->pair->vflag_either = force->pair->vflag_global = 1;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->nbondlist;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (inum > 0) {
      if (evflag) {
        if (eflag) {
          if (force->newton_bond) eval<1,1,1>(ifrom, ito, thr);
          else eval<1,1,0>(ifrom, ito, thr);
        } else {
          if (force->newton_bond) eval<1,0,1>(ifrom, ito, thr);
          else eval<1,0,0>(ifrom, ito, thr);
        }
      } else {
        if (force->newton_bond) eval<0,0,1>(ifrom, ito, thr);
        else eval<0,0,0>(ifrom, ito, thr);
      }
    }
    thr->timer(Timer::BOND);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void BondQuarticOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,n,m,type,itype,jtype;
  double delx,dely,delz,ebond,fbond,evdwl,fpair;
  double r,rsq,dr,r2,ra,rb,sr2,sr6;

  ebond = evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  int * const * const bondlist = neighbor->bondlist;
  const double * const * const cutsq = force->pair->cutsq;
  const int nlocal = atom->nlocal;

  for (n = nfrom; n < nto; n++) {

    // skip bond if already broken

    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;

    // if bond breaks, set type to 0
    //   both in temporary bondlist and permanent bond_type
    // if this proc owns both atoms,
    //   negate bond_type twice if other atom stores it
    // if other proc owns 2nd atom, other proc will also break bond

    if (rsq > rc[type]*rc[type]) {
      bondlist[n][2] = 0;
      for (m = 0; m < atom->num_bond[i1]; m++)
        if (atom->bond_atom[i1][m] == atom->tag[i2])
          atom->bond_type[i1][m] = 0;
      if (i2 < atom->nlocal)
        for (m = 0; m < atom->num_bond[i2]; m++)
          if (atom->bond_atom[i2][m] == atom->tag[i1])
            atom->bond_type[i2][m] = 0;
      continue;
    }

    // quartic bond
    // 1st portion is from quartic term
    // 2nd portion is from LJ term cut at 2^(1/6) with eps = sigma = 1.0

    r = sqrt(rsq);
    dr = r - rc[type];
    r2 = dr*dr;
    ra = dr - b1[type];
    rb = dr - b2[type];
    fbond = -k[type]/r * (r2*(ra+rb) + 2.0*dr*ra*rb);

    if (rsq < TWO_1_3) {
      sr2 = 1.0/rsq;
      sr6 = sr2*sr2*sr2;
      fbond += 48.0*sr6*(sr6-0.5)/rsq;
    }

    if (EFLAG) {
      ebond = k[type]*r2*ra*rb + u0[type];
      if (rsq < TWO_1_3) ebond += 4.0*sr6*(sr6-1.0) + 1.0;
    }

    // apply force to each of 2 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
    }

    if (EVFLAG) ev_tally_thr(this,i1,i2,nlocal,NEWTON_BOND,ebond,fbond,delx,dely,delz,thr);

    // subtract out pairwise contribution from 2 atoms via pair->single()
    // required since special_bond = 1,1,1
    // tally energy/virial in pair, using newton_bond as newton flag

    itype = atom->type[i1];
    jtype = atom->type[i2];

    if (rsq < cutsq[itype][jtype]) {
      evdwl = -force->pair->single(i1,i2,itype,jtype,rsq,1.0,1.0,fpair);
      fpair = -fpair;

      if (NEWTON_BOND || i1 < nlocal) {
        f[i1][0] += delx*fpair;
        f[i1][1] += dely*fpair;
        f[i1][2] += delz*fpair;
      }
      if (NEWTON_BOND || i2 < nlocal) {
        f[i2][0] -= delx*fpair;
        f[i2][1] -= dely*fpair;
        f[i2][2] -= delz*fpair;
      }

      if (EVFLAG) ev_tally_thr(force->pair,i1,i2,nlocal,NEWTON_BOND,
                               evdwl,0.0,fpair,delx,dely,delz,thr);
    }
  }
}

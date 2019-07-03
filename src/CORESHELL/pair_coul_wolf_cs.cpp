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

#include "pair_coul_wolf_cs.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define EPSILON 1.0e-20

/* ---------------------------------------------------------------------- */

PairCoulWolfCS::PairCoulWolfCS(LAMMPS *lmp) : PairCoulWolf( lmp )
{
   single_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairCoulWolfCS::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,forcecoul,factor_coul;
  double prefactor;
  double r;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double erfcc,erfcd,v_sh,dvdrr,e_self,e_shift,f_shift,qisq;

  ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  // self and shifted coulombic energy

  e_self = v_sh = 0.0;
  e_shift = erfc(alf*cut_coul)/cut_coul;
  f_shift = -(e_shift+ 2.0*alf/MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) /
    cut_coul;

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
    jlist = firstneigh[i];
    jnum = numneigh[i];

    qisq = qtmp*qtmp;
    e_self = -(e_shift/2.0 + alf/MY_PIS) * qisq*qqrd2e;
    if (evflag) ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_coulsq) {
        rsq += EPSILON;
        // Add EPISLON for case: r = 0; Interaction must be removed
        // by special bond
        r = sqrt(rsq);
        prefactor = qqrd2e*qtmp*q[j]/r;
        erfcc = erfc(alf*r);
        erfcd = exp(-alf*alf*r*r);
        v_sh = (erfcc - e_shift*r) * prefactor;
        dvdrr = (erfcc/rsq + 2.0*alf/MY_PIS * erfcd/r) + f_shift;
        forcecoul = dvdrr*rsq*prefactor;
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        fpair = forcecoul / rsq;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          ecoul = v_sh;
          if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
        } else ecoul = 0.0;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             0.0,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* NOTES
Using erfc and expmsq provided by math_special.h

See: http://lammps.sandia.gov/threads/msg61934.html
*/

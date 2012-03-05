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

#include "bond_morse_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "domain.h"

#include <math.h>

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondMorseOMP::BondMorseOMP(class LAMMPS *lmp)
  : BondMorse(lmp), ThrOMP(lmp,THR_BOND)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void BondMorseOMP::compute(int eflag, int vflag)
{

  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = 0;

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
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

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

    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void BondMorseOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,ralpha;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const * const bondlist = neighbor->bondlist;
  const int nlocal = atom->nlocal;

  for (n = nfrom; n < nto; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);

    dr = r - r0[type];
    ralpha = exp(-alpha[type]*dr);

    // force & energy

    if (r > 0.0) fbond = -2.0*d0[type]*alpha[type]*(1-ralpha)*ralpha/r;
    else fbond = 0.0;

    if (EFLAG) ebond = d0[type]*(1-ralpha)*(1-ralpha);

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

    if (EVFLAG) ev_tally_thr(this,i1,i2,nlocal,NEWTON_BOND,
			     ebond,fbond,delx,dely,delz,thr);
  }
}

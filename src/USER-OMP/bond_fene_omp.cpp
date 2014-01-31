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

#include "bond_fene_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "domain.h"
#include "error.h"
#include "update.h"

#include <math.h>

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondFENEOMP::BondFENEOMP(class LAMMPS *lmp)
  : BondFENE(lmp), ThrOMP(lmp,THR_BOND)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void BondFENEOMP::compute(int eflag, int vflag)
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
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void BondFENEOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r0sq,rlogarg,sr2,sr6;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int3_t * _noalias const bondlist = (int3_t *) neighbor->bondlist[0];
  const int nlocal = atom->nlocal;
  const int tid = thr->get_tid();

  for (n = nfrom; n < nto; n++) {
    i1 = bondlist[n].a;
    i2 = bondlist[n].b;
    type = bondlist[n].t;

    delx = x[i1].x - x[i2].x;
    dely = x[i1].y - x[i2].y;
    delz = x[i1].z - x[i2].z;

    rsq = delx*delx + dely*dely + delz*delz;
    r0sq = r0[type] * r0[type];
    rlogarg = 1.0 - rsq/r0sq;

    // if r -> r0, then rlogarg < 0.0 which is an error
    // issue a warning and reset rlogarg = epsilon
    // if r > 2*r0 something serious is wrong, abort

    if (rlogarg < 0.1) {
      char str[128];

      sprintf(str,"FENE bond too long: " BIGINT_FORMAT " "
              TAGINT_FORMAT " " TAGINT_FORMAT " %g",
              update->ntimestep,atom->tag[i1],atom->tag[i2],sqrt(rsq));
      error->warning(FLERR,str,0);

      if (check_error_thr((rlogarg <= -3.0),tid,FLERR,"Bad FENE bond"))
        return;

      rlogarg = 0.1;
    }

    fbond = -k[type]/rlogarg;

    // force from LJ term

    if (rsq < TWO_1_3*sigma[type]*sigma[type]) {
      sr2 = sigma[type]*sigma[type]/rsq;
      sr6 = sr2*sr2*sr2;
      fbond += 48.0*epsilon[type]*sr6*(sr6-0.5)/rsq;
    }

    // energy

    if (EFLAG) {
      ebond = -0.5 * k[type]*r0sq*log(rlogarg);
      if (rsq < TWO_1_3*sigma[type]*sigma[type])
        ebond += 4.0*epsilon[type]*sr6*(sr6-1.0) + epsilon[type];
    }

    // apply force to each of 2 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1].x += delx*fbond;
      f[i1].y += dely*fbond;
      f[i1].z += delz*fbond;
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2].x -= delx*fbond;
      f[i2].y -= dely*fbond;
      f[i2].z -= delz*fbond;
    }

    if (EVFLAG) ev_tally_thr(this,i1,i2,nlocal,NEWTON_BOND,
                             ebond,fbond,delx,dely,delz,thr);
  }
}

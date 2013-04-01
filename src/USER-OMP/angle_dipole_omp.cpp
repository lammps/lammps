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

#include "angle_dipole_omp.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "neighbor.h"
#include "domain.h"

#include "math_const.h"

#include <math.h>

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleDipoleOMP::AngleDipoleOMP(class LAMMPS *lmp)
  : AngleDipole(lmp), ThrOMP(lmp,THR_ANGLE)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void AngleDipoleOMP::compute(int eflag, int vflag)
{

  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = 0;

  if (!force->newton_bond)
    error->all(FLERR,"'newton' flag for bonded interactions must be 'on'");

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->nanglelist;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (inum > 0) {
      if (eflag)
        eval<1>(ifrom, ito, thr);
      else
        eval<0>(ifrom, ito, thr);
    }
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region

}

template <int EFLAG>
void AngleDipoleOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int iRef,iDip,iDummy,n,type;
  double delx,dely,delz;
  double eangle,tangle;
  double r,cosGamma,deltaGamma,kdg,rmu;

  const double * const * const x = atom->x;   // position vector
  const double * const * const mu = atom->mu; // point-dipole components and moment magnitude
  double * const * const torque = thr->get_torque();
  const int * const * const anglelist = neighbor->anglelist;
  const int nlocal = atom->nlocal;
  const double f1[3] = {0.0, 0.0, 0.0};
  const double f3[3] = {0.0, 0.0, 0.0};

  for (n = nfrom; n < nto; n++) {
    iDip = anglelist[n][0]; // dipole whose orientation is to be restrained
    iRef = anglelist[n][1]; // reference atom toward which dipole will point
    iDummy = anglelist[n][2]; // dummy atom - irrelevant to the interaction
    type = anglelist[n][3];

    delx = x[iRef][0] - x[iDip][0];
    dely = x[iRef][1] - x[iDip][1];
    delz = x[iRef][2] - x[iDip][2];

    r = sqrt(delx*delx + dely*dely + delz*delz);

    rmu = r * mu[iDip][3];
    cosGamma = (mu[iDip][0]*delx+mu[iDip][1]*dely+mu[iDip][2]*delz) / rmu;
    deltaGamma = cosGamma - cos(gamma0[type]);
    kdg = k[type] * deltaGamma;

    if (EFLAG) eangle = kdg * deltaGamma; // energy

    tangle = 2.0 * kdg / rmu;

    torque[iDip][0] += tangle * (dely*mu[iDip][2] - delz*mu[iDip][1]);
    torque[iDip][1] += tangle * (delz*mu[iDip][0] - delx*mu[iDip][2]);
    torque[iDip][2] += tangle * (delx*mu[iDip][1] - dely*mu[iDip][0]);

    if (EFLAG) // tally energy (virial=0 because force=0)
      ev_tally_thr(this,iRef,iDip,iDummy,nlocal,/* NEWTON_BOND */ 1,
                   eangle,f1,f3,0.0,0.0,0.0,0.0,0.0,0.0,thr);
  }
}

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

#include "angle_sdk_omp.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"

#include <math.h>

#include "lj_sdk_common.h"

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace LJSDKParms;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleSDKOMP::AngleSDKOMP(class LAMMPS *lmp)
  : AngleSDK(lmp), ThrOMP(lmp,THR_ANGLE)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void AngleSDKOMP::compute(int eflag, int vflag)
{

  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = 0;

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
void AngleSDKOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2,delx3,dely3,delz3;
  double eangle,f1[3],f3[3],e13,f13;
  double dtheta,tk;
  double rsq1,rsq2,rsq3,r1,r2,c,s,a,a11,a12,a22;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int4_t * _noalias const anglelist = (int4_t *) neighbor->anglelist[0];
  const int nlocal = atom->nlocal;

  for (n = nfrom; n < nto; n++) {
    i1 = anglelist[n].a;
    i2 = anglelist[n].b;
    i3 = anglelist[n].c;
    type = anglelist[n].t;

    // 1st bond

    delx1 = x[i1].x - x[i2].x;
    dely1 = x[i1].y - x[i2].y;
    delz1 = x[i1].z - x[i2].z;

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3].x - x[i2].x;
    dely2 = x[i3].y - x[i2].y;
    delz2 = x[i3].z - x[i2].z;

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // angle (cos and sin)

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // 1-3 LJ interaction.
    // we only want to use the repulsive part,
    // and it can be scaled (or off).
    // so this has to be done here and not in the
    // general non-bonded code.

    f13 = e13 = delx3 = dely3 = delz3 = 0.0;

    if (repflag) {

      delx3 = x[i1].x - x[i3].x;
      dely3 = x[i1].y - x[i3].y;
      delz3 = x[i1].z - x[i3].z;
      rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;

      const int type1 = atom->type[i1];
      const int type3 = atom->type[i3];

      if (rsq3 < rminsq[type1][type3]) {
        const int ljt = lj_type[type1][type3];
        const double r2inv = 1.0/rsq3;

        if (ljt == LJ12_4) {
          const double r4inv=r2inv*r2inv;

          f13 = r4inv*(lj1[type1][type3]*r4inv*r4inv - lj2[type1][type3]);
          if (EFLAG) e13 = r4inv*(lj3[type1][type3]*r4inv*r4inv - lj4[type1][type3]);

        } else if (ljt == LJ9_6) {
          const double r3inv = r2inv*sqrt(r2inv);
          const double r6inv = r3inv*r3inv;

          f13 = r6inv*(lj1[type1][type3]*r3inv - lj2[type1][type3]);
          if (EFLAG) e13 = r6inv*(lj3[type1][type3]*r3inv - lj4[type1][type3]);

        } else if (ljt == LJ12_6) {
          const double r6inv = r2inv*r2inv*r2inv;

          f13 = r6inv*(lj1[type1][type3]*r6inv - lj2[type1][type3]);
          if (EFLAG) e13 = r6inv*(lj3[type1][type3]*r6inv - lj4[type1][type3]);
        }

        // make sure energy is 0.0 at the cutoff.
        if (EFLAG) e13 -= emin[type1][type3];

        f13 *= r2inv;
      }
    }

    // force & energy

    dtheta = acos(c) - theta0[type];
    tk = k[type] * dtheta;

    if (EFLAG) eangle = tk*dtheta;

    a = -2.0 * tk * s;
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of the 3 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1].x += f1[0] + f13*delx3;
      f[i1].y += f1[1] + f13*dely3;
      f[i1].z += f1[2] + f13*delz3;
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2].x -= f1[0] + f3[0];
      f[i2].y -= f1[1] + f3[1];
      f[i2].z -= f1[2] + f3[2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3].x += f3[0] - f13*delx3;
      f[i3].y += f3[1] - f13*dely3;
      f[i3].z += f3[2] - f13*delz3;
    }

    if (EVFLAG) {
      ev_tally_thr(this,i1,i2,i3,nlocal,NEWTON_BOND,eangle,f1,f3,
                   delx1,dely1,delz1,delx2,dely2,delz2,thr);
      if (repflag) ev_tally13_thr(this,i1,i3,nlocal,NEWTON_BOND,
                                  e13,f13,delx3,dely3,delz3,thr);
    }
  }
}

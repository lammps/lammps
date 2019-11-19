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

#include "angle_cosine_periodic_omp.h"
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "timer.h"
#include "math_special.h"

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathSpecial;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleCosinePeriodicOMP::AngleCosinePeriodicOMP(class LAMMPS *lmp)
  : AngleCosinePeriodic(lmp), ThrOMP(lmp,THR_ANGLE)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void AngleCosinePeriodicOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

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
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, cvatom, thr);

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
void AngleCosinePeriodicOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i,i1,i2,i3,n,m,type,b_factor;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double rsq1,rsq2,r1,r2,c,a,a11,a12,a22;
  double tn,tn_1,tn_2,un,un_1,un_2;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int4_t * _noalias const anglelist = (int4_t *) neighbor->anglelist[0];
  const int nlocal = atom->nlocal;
  eangle = 0.0;

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

    // c = cosine of angle

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    m = multiplicity[type];
    b_factor = b[type];

    // cos(n*x) = Tn(cos(x))
    // Tn(x) = Chebyshev polynomials of the first kind: T_0 = 1, T_1 = x, ...
    // recurrence relationship:
    // Tn(x) = 2*x*T[n-1](x) - T[n-2](x) where T[-1](x) = 0
    // also, dTn(x)/dx = n*U[n-1](x)
    // where Un(x) = 2*x*U[n-1](x) - U[n-2](x) and U[-1](x) = 0
    // finally need to handle special case for n = 1

    tn = 1.0;
    tn_1 = 1.0;
    tn_2 = 0.0;
    un = 1.0;
    un_1 = 2.0;
    un_2 = 0.0;

    // force & energy

    tn_2 = c;
    for (i = 1; i <= m; i++) {
      tn = 2*c*tn_1 - tn_2;
      tn_2 = tn_1;
      tn_1 = tn;
    }

    for (i = 2; i <= m; i++) {
      un = 2*c*un_1 - un_2;
      un_2 = un_1;
      un_1 = un;
    }
    tn = b_factor*powsign(m)*tn;
    un = b_factor*powsign(m)*m*un;

    if (EFLAG) eangle = 2*k[type]*(1.0 - tn);

    a = -k[type]*un;
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of 3 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1].x += f1[0];
      f[i1].y += f1[1];
      f[i1].z += f1[2];
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2].x -= f1[0] + f3[0];
      f[i2].y -= f1[1] + f3[1];
      f[i2].z -= f1[2] + f3[2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3].x += f3[0];
      f[i3].y += f3[1];
      f[i3].z += f3[2];
    }

    if (EVFLAG) ev_tally_thr(this,i1,i2,i3,nlocal,NEWTON_BOND,eangle,f1,f3,
                             delx1,dely1,delz1,delx2,dely2,delz2,thr);
  }
}

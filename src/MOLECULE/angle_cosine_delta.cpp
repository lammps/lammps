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
   Contributing author: Axel Kohlmeyer (Temple U), akohlmey at gmail.com
------------------------------------------------------------------------- */

#include "angle_cosine_delta.h"
#include <cmath>
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"

using namespace LAMMPS_NS;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleCosineDelta::AngleCosineDelta(LAMMPS *lmp) : AngleCosineSquared(lmp) {}

/* ---------------------------------------------------------------------- */

void AngleCosineDelta::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2,theta,dtheta,dcostheta,tk;
  double eangle,f1[3],f3[3];
  double rsq1,rsq2,r1,r2,c,a,cot,a11,a12,a22,b11,b12,b22,c0,s0,s;

  eangle = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // angle (cos and sin)

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    theta = acos(c);

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    cot = c/s;

    // force & energy

    dtheta = theta - theta0[type];
    dcostheta = cos(dtheta);
    tk = k[type] * (1.0-dcostheta);

    if (eflag) eangle = tk;

    a = -k[type];

    // expand dtheta for cos and sin contribution to force

    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    b11 = -a*c*cot / rsq1;
    b12 = a*cot / (r1*r2);
    b22 = -a*c*cot / rsq2;

    c0 = cos(theta0[type]);
    s0 = sin(theta0[type]);

    f1[0] = (a11*delx1 + a12*delx2)*c0 + (b11*delx1 + b12*delx2)*s0;
    f1[1] = (a11*dely1 + a12*dely2)*c0 + (b11*dely1 + b12*dely2)*s0;
    f1[2] = (a11*delz1 + a12*delz2)*c0 + (b11*delz1 + b12*delz2)*s0;
    f3[0] = (a22*delx2 + a12*delx1)*c0 + (b22*delx2 + b12*delx1)*s0;
    f3[1] = (a22*dely2 + a12*dely1)*c0 + (b22*dely2 + b12*dely1)*s0;
    f3[2] = (a22*delz2 + a12*delz1)*c0 + (b22*delz2 + b12*delz1)*s0;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

double AngleCosineDelta::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double theta = acos(c);
  double dtheta = theta - theta0[type];
  double dcostheta = cos(dtheta);
  double tk = k[type] * (1.0-dcostheta);
  return tk;
}

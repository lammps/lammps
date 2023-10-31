/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Federica Lodesani
   (based on nb3b harmonic by Todd R. Zeitler and Stillinger-Weber pair style)
------------------------------------------------------------------------- */

#include "pair_nb3b_screened.h"

#include <cmath>

#define SMALL 0.001

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairNb3bScreened::PairNb3bScreened(LAMMPS *lmp) : PairNb3bHarmonic(lmp)
{
  variant = SCREENED;
}

/* ---------------------------------------------------------------------- */

void PairNb3bScreened::threebody(Param *paramij, Param *paramik, Param *paramijk, double rsq1,
                                 double rsq2, double *delr1, double *delr2, double *fj, double *fk,
                                 int eflag, double &eng)
{
  double dtheta, tk;
  double r1, r2, c, s, a, a11, a12, a22;
  double scr, t00, rho1inv, rho2inv;
  double ratio1, ratio2;

  // angle (cos and sin)

  r1 = sqrt(rsq1);
  r2 = sqrt(rsq2);

  c = delr1[0] * delr2[0] + delr1[1] * delr2[1] + delr1[2] * delr2[2];
  c /= r1 * r2;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  s = sqrt(1.0 - c * c);
  if (s < SMALL) s = SMALL;
  s = 1.0 / s;

  // force & energy

  // Harmonic function multiplied by a screening function
  //
  // Uijk=k/2(theta-theta0)**2 * exp[-(rij/rhoij+rik/rhoik)]
  //
  rho1inv = paramij->invrho;
  rho2inv = paramik->invrho;
  scr = exp(-r1 * rho1inv - r2 * rho2inv);

  dtheta = acos(c) - paramijk->theta0;
  tk = paramijk->k_theta * dtheta * scr;
  t00 = tk * dtheta;

  if (eflag) eng = t00;

  a = -2.0 * tk * s;
  a11 = a * c / rsq1;
  a12 = -a / (r1 * r2);
  a22 = a * c / rsq2;
  ratio1 = rho1inv / r1;
  ratio2 = rho2inv / r2;

  fj[0] = a11 * delr1[0] + a12 * delr2[0] + t00 * ratio1 * delr1[0];
  fj[1] = a11 * delr1[1] + a12 * delr2[1] + t00 * ratio1 * delr1[1];
  fj[2] = a11 * delr1[2] + a12 * delr2[2] + t00 * ratio1 * delr1[2];
  fk[0] = a22 * delr2[0] + a12 * delr1[0] + t00 * ratio2 * delr2[0];
  fk[1] = a22 * delr2[1] + a12 * delr1[1] + t00 * ratio2 * delr2[1];
  fk[2] = a22 * delr2[2] + a12 * delr1[2] + t00 * ratio2 * delr2[2];
}

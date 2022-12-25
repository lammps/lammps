// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ----------------------------------------------------------------------- */

#include "manifold_spine.h"

#include <cmath>
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace user_manifold;



manifold_spine::manifold_spine( LAMMPS *lmp, int /*argc*/, char **/*argv*/ )
  : manifold(lmp)
{
  power = 4;
}



manifold_spine_two::manifold_spine_two( LAMMPS *lmp, int narg, char **argv )
  : manifold_spine( lmp, narg, argv )
{
  power = 2;
}



/*
 * Equation for spine is:
 *
 * -x^2 - y^2 + (a^2 - (z/c)^2)*( 1 + (A*sin(const(z) *z^2))^4 )
 * const(z) = (z < 0) ? B2 : B
 * params[5] = { a, A, B, B2, c }
 */

double manifold_spine::g_and_n( const double *x, double *nn )
{
  double a  = params[0];
  double A  = params[1];
  double B  = params[2];
  double B2 = params[3];
  double c  = params[4];

  double cc, BB, z2, xy2;
  double c2, As, Ac, azc, Apart;
  double AMs, AMc;
  double dazc, dAMs;

  if (x[2] > 0) {
    BB = B;
    cc = c;
  } else {
    BB = B2;
    cc = 1.0;
  }

  xy2 = x[0]*x[0] + x[1]*x[1];
  z2  = x[2]*x[2];
  c2  = 1.0/( cc*cc );

  azc = a*a - z2*c2;
  As  = sin( BB*z2 );
  Ac  = cos( BB*z2 );
  As *= A;
  Ac *= A;

  Apart = MathSpecial::powint( As, power-1 ); // Much more efficient! =D
  //Apart = pow( As, power-1 );
  AMs   = Apart*As;
  AMc   = Apart*Ac;

  dAMs = power * AMc * 2.0*BB*x[2];
  dazc = -2.0*x[2]*c2;

  nn[0] = -2*x[0];
  nn[1] = -2*x[1];
  nn[2] = azc * dAMs + ( 1.0 + AMs ) * dazc;

  return -xy2 + azc * ( 1.0 + AMs );
}




void manifold_spine::n( const double *x, double *nn )
{
  double a  = params[0];
  double A  = params[1];
  double B  = params[2];
  double B2 = params[3];
  double c  = params[4];

  double cc, BB, z2;
  double c2, As, Ac, azc, Apart;
  double AMs, AMc;
  double dazc, dAMs;

  if (x[2] > 0) {
    BB = B;
    cc = c;
  } else {
    BB = B2;
    cc = 1.0;
  }

  z2  = x[2]*x[2];
  c2  = 1.0/( cc*cc );

  azc = a*a - z2*c2;
  As  = sin( BB*z2 );
  Ac  = cos( BB*z2 );
  As *= A;
  Ac *= A;

  Apart = MathSpecial::powint( As, power-1 ); // Much more efficient! =D
  //Apart = pow( As, power-1 );
  AMs   = Apart*As;
  AMc   = Apart*Ac;

  dAMs = power * AMc * 2.0*BB*x[2];
  dazc = -2.0*x[2]*c2;

  nn[0] = -2*x[0];
  nn[1] = -2*x[1];
  nn[2] = azc * dAMs + ( 1.0 + AMs ) * dazc;

}


double manifold_spine::g( const double *x )
{
  double a  = params[0];
  double A  = params[1];
  double B  = params[2];
  double B2 = params[3];
  double c  = params[4];

  double cc, BB, z2, xy2;
  double c2, As, azc, Apart;
  double AMs;

  if (x[2] > 0) {
    BB = B;
    cc = c;
  } else {
    BB = B2;
    cc = 1.0;
  }

  xy2 = x[0]*x[0] + x[1]*x[1];
  z2  = x[2]*x[2];
  c2  = 1.0/( cc*cc );

  azc = a*a - z2*c2;
  As  = sin( BB*z2 );
  As *= A;

  Apart = MathSpecial::powint( As, power-1 );
  AMs   = Apart*As;
  return -xy2 + azc * ( 1.0 + AMs );

}



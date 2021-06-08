// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ----------------------------------------------------------------------- */

#include <cmath>
#include "manifold_torus.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace user_manifold;


manifold_torus::manifold_torus( LAMMPS *lmp, int /*argc*/, char **/*argv*/ ) : manifold(lmp)
{}


double manifold_torus::g( const double *x )
{
  double R = params[0];
  double r = params[1];
  if (R < r) {
    error->all(FLERR,"Large radius < small radius!");
  }

  double rad = sqrt(x[0]*x[0] + x[1]*x[1]);
  double c = R - rad;
  return c*c + x[2]*x[2] - r*r;
}

void   manifold_torus::n( const double *x, double *n )
{
  double R = params[0];
  double r = params[1];
  if (R < r) {
    error->all(FLERR,"Large radius < small radius!");
  }

  double rad = sqrt(x[0]*x[0] + x[1]*x[1]);
  double c = R - rad;
  double fac = c / rad;
  n[0] = -2.0 * fac * x[0];
  n[1] = -2.0 * fac * x[1];
  n[2] = 2.0 * x[2];
}

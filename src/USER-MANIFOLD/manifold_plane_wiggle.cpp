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

#include "manifold_plane_wiggle.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace user_manifold;

manifold_plane_wiggle::manifold_plane_wiggle( LAMMPS *lmp, int /*argc*/, char **/*argv*/ ) :
  manifold(lmp)
{}


double manifold_plane_wiggle::g( const double *x )
{
  double a = params[0];
  double w = params[1];
  return x[2] - a*sin(w*x[0]);
}


void manifold_plane_wiggle::n( const double *x, double *n )
{
  double a = params[0];
  double w = params[1];
  n[2] = 1;
  n[1] = 0.0;
  n[0] = -a*w*cos(w*x[0]);
}

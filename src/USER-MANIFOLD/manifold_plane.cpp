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

#include "manifold_plane.h"

using namespace LAMMPS_NS;

using namespace user_manifold;

manifold_plane::manifold_plane( LAMMPS *lmp, int /*argc*/, char **/*argv*/ ) :
  manifold(lmp)
{}

double manifold_plane::g( const double *x )
{
  double a  = params[0], b  = params[1], c  = params[2];
  double x0 = params[3], y0 = params[4], z0 = params[5];
  return a*(x[0] - x0) + b*(x[1] - y0) + c*(x[2] - z0);
}


void manifold_plane::n( const double * /*x*/, double *n )
{
  n[0] = params[0];
  n[1] = params[1];
  n[2] = params[2];
}

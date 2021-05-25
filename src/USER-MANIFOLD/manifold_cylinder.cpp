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

#include "manifold_cylinder.h"

using namespace LAMMPS_NS;

using namespace user_manifold;


manifold_cylinder::manifold_cylinder( LAMMPS *lmp, int /*argc*/,
                                      char **/*argv*/ ) : manifold(lmp)
{}


double manifold_cylinder::g( const double *x )
{
  double R = params[0];
  double r2 = x[0]*x[0] + x[1]*x[1];
  return R*R - r2;
}

void manifold_cylinder::n( const double *x, double *n )
{
  n[0] = -2*x[0];
  n[1] = -2*x[1];
  n[2] = 0.0;
}



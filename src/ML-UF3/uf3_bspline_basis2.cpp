// clang-format off
/* ----------------------------------------------------------------------
   lammps - large-scale atomic/molecular massively parallel simulator
   https://www.lammps.org/, sandia national laboratories
   lammps development team: developers@lammps.org

   copyright (2003) sandia corporation.  under the terms of contract
   de-ac04-94al85000 with sandia corporation, the u.s. government retains
   certain rights in this software.  this software is distributed under
   the gnu general public license.

   see the readme file in the top-level lammps directory.
------------------------------------------------------------------------- */

#include "uf3_bspline_basis2.h"

#include "math_special.h"

using namespace LAMMPS_NS;
using MathSpecial::square;

// Constructor
// Initializes coefficients and knots
// Requires [knots] to have length 4
uf3_bspline_basis2::uf3_bspline_basis2(LAMMPS *ulmp, const double *knots, double coefficient)
{
  lmp = ulmp;

  double c0, c1, c2;

  c0 = coefficient
    * (square(knots[0])
       / (square(knots[0]) - knots[0] * knots[1] - knots[0] * knots[2] + knots[1] * knots[2]));
  c1 = coefficient *
      (-2.0 * knots[0] /
       (square(knots[0]) - knots[0] * knots[1] - knots[0] * knots[2] + knots[1] * knots[2]));
  c2 = coefficient *
      (1.0 / (square(knots[0]) - knots[0] * knots[1] - knots[0] * knots[2] + knots[1] * knots[2]));
  //constants.push_back(c0);
  //constants.push_back(c1);
  //constants.push_back(c2);
  constants[0] = c0;
  constants[1] = c1;
  constants[2] = c2;
  c0 = coefficient *
      (-knots[1] * knots[3] /
           (square(knots[1]) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) -
       knots[0] * knots[2] /
           (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + square(knots[2])));
  c1 = coefficient *
      (knots[1] /
           (square(knots[1]) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) +
       knots[3] /
           (square(knots[1]) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) +
       knots[0] /
           (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + square(knots[2])) +
       knots[2] /
           (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + square(knots[2])));
  c2 = coefficient *
      (-1.0 / (square(knots[1]) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) -
       1.0 / (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + square(knots[2])));
  //constants.push_back(c0);
  //constants.push_back(c1);
  //constants.push_back(c2);
  constants[3] = c0;
  constants[4] = c1;
  constants[5] = c2;
  c0 = coefficient *
      (square(knots[3]) /
       (knots[1] * knots[2] - knots[1] * knots[3] - knots[2] * knots[3] + square(knots[3])));
  c1 = coefficient *
      (-2.0 * knots[3] /
       (knots[1] * knots[2] - knots[1] * knots[3] - knots[2] * knots[3] + square(knots[3])));
  c2 = coefficient *
      (1.0 / (knots[1] * knots[2] - knots[1] * knots[3] - knots[2] * knots[3] + square(knots[3])));
  //constants.push_back(c0);
  //constants.push_back(c1);
  //constants.push_back(c2);
  constants[6] = c0;
  constants[7] = c1;
  constants[8] = c2;
}

uf3_bspline_basis2::~uf3_bspline_basis2() {}

// Evaluate outer-left part of spline
double uf3_bspline_basis2::eval0(double rsq, double r)
{
  return rsq * constants[2] + r * constants[1] + constants[0];
}

// Evaluate center-left part of spline
double uf3_bspline_basis2::eval1(double rsq, double r)
{
  return rsq * constants[5] + r * constants[4] + constants[3];
}

// Evaluate center-right part of spline
double uf3_bspline_basis2::eval2(double rsq, double r)
{
  return rsq * constants[8] + r * constants[7] + constants[6];
}

double uf3_bspline_basis2::memory_usage()
{
  double bytes = 0;

  bytes += (double)9*sizeof(double);

  return bytes;
}

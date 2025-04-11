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

#include "uf3_bspline_basis3.h"

#include "math_special.h"

using namespace LAMMPS_NS;
using MathSpecial::cube;
using MathSpecial::square;

// Constructor
// Initializes coefficients and knots
// [knots] needs to have length 4
uf3_bspline_basis3::uf3_bspline_basis3(LAMMPS *ulmp, const double *knots, double coefficient)
{
  lmp = ulmp;

  double c0, c1, c2, c3;

  c0 = coefficient *
      (-cube(knots[0]) /
       (-cube(knots[0]) + square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
        square(knots[0]) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  c1 = coefficient *
      (3.0 * square(knots[0]) /
       (-cube(knots[0]) + square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
        square(knots[0]) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  c2 = coefficient *
      (-3.0 * knots[0] /
       (-cube(knots[0]) + square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
        square(knots[0]) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  c3 = coefficient *
      (1.0 /
       (-cube(knots[0]) + square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
        square(knots[0]) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  //constants.push_back(c0);
  //constants.push_back(c1);
  //constants.push_back(c2);
  //constants.push_back(c3);
  constants[0] = c0;
  constants[1] = c1;
  constants[2] = c2;
  constants[3] = c3;
  c0 = coefficient *
      (square(knots[1]) * knots[4] /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) +
       square(knots[0]) * knots[2] /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) +
       knots[0] * knots[1] * knots[3] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])));
  c1 = coefficient *
      (-square(knots[1]) /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) -
       2.0 * knots[1] * knots[4] /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) -
       square(knots[0]) /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) -
       2.0 * knots[0] * knots[2] /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) -
       knots[0] * knots[1] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])) -
       knots[0] * knots[3] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])) -
       knots[1] * knots[3] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])));
  c2 = coefficient *
      (2.0 * knots[1] /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) +
       knots[4] /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) +
       2.0 * knots[0] /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) +
       knots[2] /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) +
       knots[0] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])) +
       knots[1] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])) +
       knots[3] /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])));
  c3 = coefficient *
      (-1.0 /
           (-cube(knots[1]) + square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            square(knots[1]) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) -
       1.0 /
           (-square(knots[0]) * knots[1] + square(knots[0]) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * square(knots[2]) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + square(knots[2]) * knots[3]) -
       1.0 /
           (-knots[0] * square(knots[1]) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            square(knots[1]) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * square(knots[3]) + knots[2] * square(knots[3])));
  //constants.push_back(c0);
  //constants.push_back(c1);
  //constants.push_back(c2);
  //constants.push_back(c3);
  constants[4] = c0;
  constants[5] = c1;
  constants[6] = c2;
  constants[7] = c3;
  c0 = coefficient *
      (-knots[0] * square(knots[3]) /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) -
       knots[1] * knots[3] * knots[4] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) -
       knots[2] * square(knots[4]) /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])));
  c1 = coefficient *
      (2.0 * knots[0] * knots[3] /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) +
       square(knots[3]) /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) +
       knots[1] * knots[3] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) +
       knots[1] * knots[4] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) +
       knots[3] * knots[4] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) +
       2.0 * knots[2] * knots[4] /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])) +
       square(knots[4]) /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])));
  c2 = coefficient *
      (-knots[0] /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) -
       2.0 * knots[3] /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) -
       knots[1] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) -
       knots[3] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) -
       knots[4] /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) -
       knots[2] /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])) -
       2.0 * knots[4] /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])));
  c3 = coefficient *
      (1.0 /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * square(knots[3]) +
            knots[1] * knots[2] * knots[3] - knots[1] * square(knots[3]) -
            knots[2] * square(knots[3]) + cube(knots[3])) +
       1.0 /
           (-square(knots[1]) * knots[2] + square(knots[1]) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * square(knots[3]) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + square(knots[3]) * knots[4]) +
       1.0 /
           (-knots[1] * square(knots[2]) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            square(knots[2]) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * square(knots[4]) + knots[3] * square(knots[4])));
  //constants.push_back(c0);
  //constants.push_back(c1);
  //constants.push_back(c2);
  //constants.push_back(c3);
  constants[8] = c0;
  constants[9] = c1;
  constants[10] = c2;
  constants[11] = c3;
  c0 = coefficient *
      (cube(knots[4]) /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * square(knots[4]) +
        knots[2] * knots[3] * knots[4] - knots[2] * square(knots[4]) - knots[3] * square(knots[4]) +
        cube(knots[4])));
  c1 = coefficient *
      (-3.0 * square(knots[4]) /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * square(knots[4]) +
        knots[2] * knots[3] * knots[4] - knots[2] * square(knots[4]) - knots[3] * square(knots[4]) +
        cube(knots[4])));
  c2 = coefficient *
      (3.0 * knots[4] /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * square(knots[4]) +
        knots[2] * knots[3] * knots[4] - knots[2] * square(knots[4]) - knots[3] * square(knots[4]) +
        cube(knots[4])));
  c3 = coefficient *
      (-1.0 /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * square(knots[4]) +
        knots[2] * knots[3] * knots[4] - knots[2] * square(knots[4]) - knots[3] * square(knots[4]) +
        cube(knots[4])));
  //constants.push_back(c0);
  //constants.push_back(c1);
  //constants.push_back(c2);
  //constants.push_back(c3);
  constants[12] = c0;
  constants[13] = c1;
  constants[14] = c2;
  constants[15] = c3;
}

uf3_bspline_basis3::~uf3_bspline_basis3() {}

// Evaluate outer-left part of spline
double uf3_bspline_basis3::eval0(double rth, double rsq, double r)
{
  return rth * constants[3] + rsq * constants[2] + r * constants[1] + constants[0];
}

// Evaluate center-left part of spline
double uf3_bspline_basis3::eval1(double rth, double rsq, double r)
{
  return rth * constants[7] + rsq * constants[6] + r * constants[5] + constants[4];
}

// Evaluate center-right part of spline
double uf3_bspline_basis3::eval2(double rth, double rsq, double r)
{
  return rth * constants[11] + rsq * constants[10] + r * constants[9] + constants[8];
}

// Evaluate outer-right part of spline
double uf3_bspline_basis3::eval3(double rth, double rsq, double r)
{
  return rth * constants[15] + rsq * constants[14] + r * constants[13] + constants[12];
}

double uf3_bspline_basis3::memory_usage()
{
  double bytes = 0;

  bytes += (double)16*sizeof(double);

  return bytes;
}

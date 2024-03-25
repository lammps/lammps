#include "uf3_bspline_basis3.h"

#include "utils.h"
#include <vector>

using namespace LAMMPS_NS;

// Constructor
// Initializes coefficients and knots
// [knots] needs to have length 4
uf3_bspline_basis3::uf3_bspline_basis3(LAMMPS *ulmp, const double *knots, double coefficient)
{
  lmp = ulmp;

  double c0, c1, c2, c3;

  c0 = coefficient *
      (-pow(knots[0], 3) /
       (-pow(knots[0], 3) + pow(knots[0], 2) * knots[1] + pow(knots[0], 2) * knots[2] +
        pow(knots[0], 2) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  c1 = coefficient *
      (3 * pow(knots[0], 2) /
       (-pow(knots[0], 3) + pow(knots[0], 2) * knots[1] + pow(knots[0], 2) * knots[2] +
        pow(knots[0], 2) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  c2 = coefficient *
      (-3 * knots[0] /
       (-pow(knots[0], 3) + pow(knots[0], 2) * knots[1] + pow(knots[0], 2) * knots[2] +
        pow(knots[0], 2) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  c3 = coefficient *
      (1 /
       (-pow(knots[0], 3) + pow(knots[0], 2) * knots[1] + pow(knots[0], 2) * knots[2] +
        pow(knots[0], 2) * knots[3] - knots[0] * knots[1] * knots[2] -
        knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
        knots[1] * knots[2] * knots[3]));
  constants.push_back(c0);
  constants.push_back(c1);
  constants.push_back(c2);
  constants.push_back(c3);
  c0 = coefficient *
      (pow(knots[1], 2) * knots[4] /
           (-pow(knots[1], 3) + pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            pow(knots[1], 2) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) +
       pow(knots[0], 2) * knots[2] /
           (-pow(knots[0], 2) * knots[1] + pow(knots[0], 2) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * pow(knots[2], 2) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + pow(knots[2], 2) * knots[3]) +
       knots[0] * knots[1] * knots[3] /
           (-knots[0] * pow(knots[1], 2) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            pow(knots[1], 2) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * pow(knots[3], 2) + knots[2] * pow(knots[3], 2)));
  c1 = coefficient *
      (-pow(knots[1], 2) /
           (-pow(knots[1], 3) + pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            pow(knots[1], 2) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) -
       2 * knots[1] * knots[4] /
           (-pow(knots[1], 3) + pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            pow(knots[1], 2) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) -
       pow(knots[0], 2) /
           (-pow(knots[0], 2) * knots[1] + pow(knots[0], 2) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * pow(knots[2], 2) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + pow(knots[2], 2) * knots[3]) -
       2 * knots[0] * knots[2] /
           (-pow(knots[0], 2) * knots[1] + pow(knots[0], 2) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * pow(knots[2], 2) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + pow(knots[2], 2) * knots[3]) -
       knots[0] * knots[1] /
           (-knots[0] * pow(knots[1], 2) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            pow(knots[1], 2) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * pow(knots[3], 2) + knots[2] * pow(knots[3], 2)) -
       knots[0] * knots[3] /
           (-knots[0] * pow(knots[1], 2) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            pow(knots[1], 2) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * pow(knots[3], 2) + knots[2] * pow(knots[3], 2)) -
       knots[1] * knots[3] /
           (-knots[0] * pow(knots[1], 2) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            pow(knots[1], 2) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * pow(knots[3], 2) + knots[2] * pow(knots[3], 2)));
  c2 = coefficient *
      (2 * knots[1] /
           (-pow(knots[1], 3) + pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            pow(knots[1], 2) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) +
       knots[4] /
           (-pow(knots[1], 3) + pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            pow(knots[1], 2) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) +
       2 * knots[0] /
           (-pow(knots[0], 2) * knots[1] + pow(knots[0], 2) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * pow(knots[2], 2) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + pow(knots[2], 2) * knots[3]) +
       knots[2] /
           (-pow(knots[0], 2) * knots[1] + pow(knots[0], 2) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * pow(knots[2], 2) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + pow(knots[2], 2) * knots[3]) +
       knots[0] /
           (-knots[0] * pow(knots[1], 2) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            pow(knots[1], 2) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * pow(knots[3], 2) + knots[2] * pow(knots[3], 2)) +
       knots[1] /
           (-knots[0] * pow(knots[1], 2) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            pow(knots[1], 2) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * pow(knots[3], 2) + knots[2] * pow(knots[3], 2)) +
       knots[3] /
           (-knots[0] * pow(knots[1], 2) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            pow(knots[1], 2) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * pow(knots[3], 2) + knots[2] * pow(knots[3], 2)));
  c3 = coefficient *
      (-1 /
           (-pow(knots[1], 3) + pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            pow(knots[1], 2) * knots[4] - knots[1] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            knots[2] * knots[3] * knots[4]) -
       1 /
           (-pow(knots[0], 2) * knots[1] + pow(knots[0], 2) * knots[2] +
            knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] -
            knots[0] * pow(knots[2], 2) - knots[0] * knots[2] * knots[3] -
            knots[1] * knots[2] * knots[3] + pow(knots[2], 2) * knots[3]) -
       1 /
           (-knots[0] * pow(knots[1], 2) + knots[0] * knots[1] * knots[2] +
            knots[0] * knots[1] * knots[3] - knots[0] * knots[2] * knots[3] +
            pow(knots[1], 2) * knots[3] - knots[1] * knots[2] * knots[3] -
            knots[1] * pow(knots[3], 2) + knots[2] * pow(knots[3], 2)));
  constants.push_back(c0);
  constants.push_back(c1);
  constants.push_back(c2);
  constants.push_back(c3);
  c0 = coefficient *
      (-knots[0] * pow(knots[3], 2) /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * pow(knots[3], 2) +
            knots[1] * knots[2] * knots[3] - knots[1] * pow(knots[3], 2) -
            knots[2] * pow(knots[3], 2) + pow(knots[3], 3)) -
       knots[1] * knots[3] * knots[4] /
           (-pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * pow(knots[3], 2) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + pow(knots[3], 2) * knots[4]) -
       knots[2] * pow(knots[4], 2) /
           (-knots[1] * pow(knots[2], 2) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            pow(knots[2], 2) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * pow(knots[4], 2) + knots[3] * pow(knots[4], 2)));
  c1 = coefficient *
      (2 * knots[0] * knots[3] /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * pow(knots[3], 2) +
            knots[1] * knots[2] * knots[3] - knots[1] * pow(knots[3], 2) -
            knots[2] * pow(knots[3], 2) + pow(knots[3], 3)) +
       pow(knots[3], 2) /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * pow(knots[3], 2) +
            knots[1] * knots[2] * knots[3] - knots[1] * pow(knots[3], 2) -
            knots[2] * pow(knots[3], 2) + pow(knots[3], 3)) +
       knots[1] * knots[3] /
           (-pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * pow(knots[3], 2) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + pow(knots[3], 2) * knots[4]) +
       knots[1] * knots[4] /
           (-pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * pow(knots[3], 2) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + pow(knots[3], 2) * knots[4]) +
       knots[3] * knots[4] /
           (-pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * pow(knots[3], 2) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + pow(knots[3], 2) * knots[4]) +
       2 * knots[2] * knots[4] /
           (-knots[1] * pow(knots[2], 2) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            pow(knots[2], 2) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * pow(knots[4], 2) + knots[3] * pow(knots[4], 2)) +
       pow(knots[4], 2) /
           (-knots[1] * pow(knots[2], 2) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            pow(knots[2], 2) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * pow(knots[4], 2) + knots[3] * pow(knots[4], 2)));
  c2 = coefficient *
      (-knots[0] /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * pow(knots[3], 2) +
            knots[1] * knots[2] * knots[3] - knots[1] * pow(knots[3], 2) -
            knots[2] * pow(knots[3], 2) + pow(knots[3], 3)) -
       2 * knots[3] /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * pow(knots[3], 2) +
            knots[1] * knots[2] * knots[3] - knots[1] * pow(knots[3], 2) -
            knots[2] * pow(knots[3], 2) + pow(knots[3], 3)) -
       knots[1] /
           (-pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * pow(knots[3], 2) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + pow(knots[3], 2) * knots[4]) -
       knots[3] /
           (-pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * pow(knots[3], 2) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + pow(knots[3], 2) * knots[4]) -
       knots[4] /
           (-pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * pow(knots[3], 2) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + pow(knots[3], 2) * knots[4]) -
       knots[2] /
           (-knots[1] * pow(knots[2], 2) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            pow(knots[2], 2) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * pow(knots[4], 2) + knots[3] * pow(knots[4], 2)) -
       2 * knots[4] /
           (-knots[1] * pow(knots[2], 2) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            pow(knots[2], 2) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * pow(knots[4], 2) + knots[3] * pow(knots[4], 2)));
  c3 = coefficient *
      (1 /
           (-knots[0] * knots[1] * knots[2] + knots[0] * knots[1] * knots[3] +
            knots[0] * knots[2] * knots[3] - knots[0] * pow(knots[3], 2) +
            knots[1] * knots[2] * knots[3] - knots[1] * pow(knots[3], 2) -
            knots[2] * pow(knots[3], 2) + pow(knots[3], 3)) +
       1 /
           (-pow(knots[1], 2) * knots[2] + pow(knots[1], 2) * knots[3] +
            knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] -
            knots[1] * pow(knots[3], 2) - knots[1] * knots[3] * knots[4] -
            knots[2] * knots[3] * knots[4] + pow(knots[3], 2) * knots[4]) +
       1 /
           (-knots[1] * pow(knots[2], 2) + knots[1] * knots[2] * knots[3] +
            knots[1] * knots[2] * knots[4] - knots[1] * knots[3] * knots[4] +
            pow(knots[2], 2) * knots[4] - knots[2] * knots[3] * knots[4] -
            knots[2] * pow(knots[4], 2) + knots[3] * pow(knots[4], 2)));
  constants.push_back(c0);
  constants.push_back(c1);
  constants.push_back(c2);
  constants.push_back(c3);
  c0 = coefficient *
      (pow(knots[4], 3) /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * pow(knots[4], 2) +
        knots[2] * knots[3] * knots[4] - knots[2] * pow(knots[4], 2) - knots[3] * pow(knots[4], 2) +
        pow(knots[4], 3)));
  c1 = coefficient *
      (-3 * pow(knots[4], 2) /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * pow(knots[4], 2) +
        knots[2] * knots[3] * knots[4] - knots[2] * pow(knots[4], 2) - knots[3] * pow(knots[4], 2) +
        pow(knots[4], 3)));
  c2 = coefficient *
      (3 * knots[4] /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * pow(knots[4], 2) +
        knots[2] * knots[3] * knots[4] - knots[2] * pow(knots[4], 2) - knots[3] * pow(knots[4], 2) +
        pow(knots[4], 3)));
  c3 = coefficient *
      (-1 /
       (-knots[1] * knots[2] * knots[3] + knots[1] * knots[2] * knots[4] +
        knots[1] * knots[3] * knots[4] - knots[1] * pow(knots[4], 2) +
        knots[2] * knots[3] * knots[4] - knots[2] * pow(knots[4], 2) - knots[3] * pow(knots[4], 2) +
        pow(knots[4], 3)));
  constants.push_back(c0);
  constants.push_back(c1);
  constants.push_back(c2);
  constants.push_back(c3);
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

  bytes += (double)constants.size()*sizeof(double);

  return bytes;
}

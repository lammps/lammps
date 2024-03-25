#include "uf3_bspline_basis2.h"

#include "utils.h"
#include <vector>

using namespace LAMMPS_NS;

// Constructor
// Initializes coefficients and knots
// Requires [knots] to have length 4
uf3_bspline_basis2::uf3_bspline_basis2(LAMMPS *ulmp, const double *knots, double coefficient)
{
  lmp = ulmp;

  double c0, c1, c2;

  c0 = coefficient *
      (pow(knots[0], 2) /
       (pow(knots[0], 2) - knots[0] * knots[1] - knots[0] * knots[2] + knots[1] * knots[2]));
  c1 = coefficient *
      (-2 * knots[0] /
       (pow(knots[0], 2) - knots[0] * knots[1] - knots[0] * knots[2] + knots[1] * knots[2]));
  c2 = coefficient *
      (1 / (pow(knots[0], 2) - knots[0] * knots[1] - knots[0] * knots[2] + knots[1] * knots[2]));
  constants.push_back(c0);
  constants.push_back(c1);
  constants.push_back(c2);
  c0 = coefficient *
      (-knots[1] * knots[3] /
           (pow(knots[1], 2) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) -
       knots[0] * knots[2] /
           (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + pow(knots[2], 2)));
  c1 = coefficient *
      (knots[1] /
           (pow(knots[1], 2) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) +
       knots[3] /
           (pow(knots[1], 2) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) +
       knots[0] /
           (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + pow(knots[2], 2)) +
       knots[2] /
           (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + pow(knots[2], 2)));
  c2 = coefficient *
      (-1 / (pow(knots[1], 2) - knots[1] * knots[2] - knots[1] * knots[3] + knots[2] * knots[3]) -
       1 / (knots[0] * knots[1] - knots[0] * knots[2] - knots[1] * knots[2] + pow(knots[2], 2)));
  constants.push_back(c0);
  constants.push_back(c1);
  constants.push_back(c2);
  c0 = coefficient *
      (pow(knots[3], 2) /
       (knots[1] * knots[2] - knots[1] * knots[3] - knots[2] * knots[3] + pow(knots[3], 2)));
  c1 = coefficient *
      (-2 * knots[3] /
       (knots[1] * knots[2] - knots[1] * knots[3] - knots[2] * knots[3] + pow(knots[3], 2)));
  c2 = coefficient *
      (1 / (knots[1] * knots[2] - knots[1] * knots[3] - knots[2] * knots[3] + pow(knots[3], 2)));
  constants.push_back(c0);
  constants.push_back(c1);
  constants.push_back(c2);
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

  bytes += (double)constants.size()*sizeof(double);

  return bytes;
}

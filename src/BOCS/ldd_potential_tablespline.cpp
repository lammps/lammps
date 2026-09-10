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
/* ------------------------------------------------------
    This file is part of the BOCS package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */
#include "ldd_potential_tablespline.h"

#include "error.h"
#include "memory.h"
#include "utils.h"

#include <cmath>

using namespace LAMMPS_NS;

namespace {
constexpr double BIG = 0.99e30;
constexpr double BIGGER = 0.999e30;
constexpr double SMALL = 1.0e-8;
}    // namespace

LddPotentialTableSpline::LddPotentialTableSpline(class LAMMPS *lmp) : LddPotential(lmp)
{
  n_coeffs = 1;
}

LddPotentialTableSpline::~LddPotentialTableSpline()
{
  if (allocated == 1) {
    memory->destroy(coeffs);
    memory->destroy(potl_table.r);
    memory->destroy(potl_table.u);
    memory->destroy(potl_table.f);
    memory->destroy(potl_table.u2);
    memory->destroy(potl_table.f2);
  }
  allocated = 0;
}

void LddPotentialTableSpline::allocate()
{
  memory->create(coeffs, n_coeffs, "ldd_potential:coeffs");
  allocated = 1;
}

/*******************************
Calculates second derivatives at knot points. stores them in y2[]. pass in
knot points x[] and y[]. pass in the number of knot points n. pass in first
derivative you want for the end points in yp1 and ypn. if yp1 and/or ypn
are >= 0.999e30, then second derivative is set to zero for that boundary.
Copied from "Numerical Recipes in C" second edition.
*******************************/
void LddPotentialTableSpline::spline(double x[], double y[], int n, double yp1, double ypn,
                                     double y2[])
{
  int i, j;
  double p, qn, sig, un;
  std::vector<double> u(n);
  if (yp1 > BIG)
    y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }

  for (i = 1; i < n - 1; i++) {
    sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }

  if (ypn > BIG)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }

  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

  for (j = n - 2; j > 0; j--) { y2[j] = y2[j] * y2[j + 1] + u[j]; }
}

void LddPotentialTableSpline::setup_potl(int ipt, int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg <= ipt + 1)
    error->all(FLERR, "Missing filename following the ldd table/spline keyword");
  read_table_file(arg[ipt + 2], true);

  spline(&(potl_table.r[0]), &(potl_table.u[0]), potl_table.n_pts, BIGGER, BIGGER,
         &(potl_table.u2[0]));
  spline(&(potl_table.r[0]), &(potl_table.f[0]), potl_table.n_pts, BIGGER, BIGGER,
         &(potl_table.f2[0]));
}

/* Evaluate cubic spline. Modified from "Numerical Recipes in C" Second Edition */
double LddPotentialTableSpline::splint(double x0, double x1, double y0, double y1, double y20,
                                       double y21, double dr, double x, double a, double b)
{
  if (fabs(x - x0) <= SMALL) return y0;
  if (fabs(x1 - x) <= SMALL) return y1;
  return (a * y0 + b * y1 + ((a * a * a - a) * y20 + (b * b * b - b) * y21) * dr * dr / 6.0);
}

double LddPotentialTableSpline::u(double rho)
{
  // Handle this case separately
  if (rho == potl_table.r[potl_table.n_pts - 1]) { return potl_table.u[potl_table.n_pts - 1]; }
  int idx = get_table_index(rho);
  double A = calc_A_table(rho, idx);

  return splint(potl_table.r[idx], potl_table.r[idx + 1], potl_table.u[idx], potl_table.u[idx + 1],
                potl_table.u2[idx], potl_table.u2[idx + 1], potl_table.dr, rho, A, 1.0 - A);
}

double LddPotentialTableSpline::f(double rho)
{
  if (rho == potl_table.r[potl_table.n_pts - 1]) { return potl_table.f[potl_table.n_pts - 1]; }
  int idx = get_table_index(rho);
  double A = calc_A_table(rho, idx);

  return splint(potl_table.r[idx], potl_table.r[idx + 1], potl_table.f[idx], potl_table.f[idx + 1],
                potl_table.f2[idx], potl_table.f2[idx + 1], potl_table.dr, rho, A, 1.0 - A);
}

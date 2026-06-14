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
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */

#include "ldd_indicator_sphere.h"

#include "math_const.h"
#include "math_special.h"
#include "memory.h"

using namespace LAMMPS_NS;
using MathConst::MY_PI;
using MathSpecial::powint;

LddIndicatorSphere::LddIndicatorSphere(class LAMMPS *lmp) : LddIndicator(lmp)
{
  n_coeffs = 4;
}

LddIndicatorSphere::~LddIndicatorSphere()
{
  if (allocated == 1) memory->destroy(coeffs);
  allocated = 0;
}

void LddIndicatorSphere::allocate()
{
  memory->create(coeffs, n_coeffs, "ldd_indicator:coeffs");
  allocated = 1;
}

void LddIndicatorSphere::init_coeffs(double a, double b, int dim)
{
  if (!allocated) allocate();
  r0 = a;
  rc = b;
  rs = (b - a) / 2.0;
  rb = (b + a) / 2.0;
  // If only C++ could use negative indices like python...
  coeffs[2] = -3.0 * (powint(rb, 4) + powint(rs, 4)) / (16.0 * powint(rs, 3)) +
      3.0 * powint(rb, 2) / (8.0 * rs);
  coeffs[0] = (powint(rb, 3) + powint(rs, 3)) / (2.0 * powint(rs, 3));
  coeffs[1] = -3.0 * (rb * rb + rs * rs) / (8.0 * powint(rs, 3));
  coeffs[3] = 1.0 / (16.0 * powint(rs, 3));
  switch (dim) {
    case 2:
      norm = 2.0 * MY_PI *
          (powint(r0, 2) / 2.0 + coeffs[2] * (rc - r0) +
           coeffs[0] / 2.0 * (powint(rc, 2) - powint(r0, 2)) +
           coeffs[1] / 3.0 * (powint(rc, 3) - powint(r0, 3)) +
           coeffs[3] / 5.0 * (powint(rc, 5) - powint(r0, 5)));
      break;
    case 3:
      norm = 4.0 * MY_PI *
          (powint(r0, 3) / 3.0 + coeffs[2] / 2.0 * (powint(rc, 2) - powint(r0, 2)) +
           coeffs[0] / 3.0 * (powint(rc, 3) - powint(r0, 3)) +
           coeffs[1] / 4.0 * (powint(rc, 4) - powint(r0, 4)) +
           coeffs[3] / 6.0 * (powint(rc, 6) - powint(r0, 6)));
      break;
  }
  invnorm = 1.0 / norm;
}

double LddIndicatorSphere::w(double r)
{
  if (r < r0) { return invnorm; }
  if (r > rc) { return 0.0; }
  return (coeffs[2] / r + coeffs[0] + coeffs[1] * r + coeffs[3] * powint(r, 3)) * invnorm;
}

double LddIndicatorSphere::wp(double r)
{
  if (r < r0) { return 0.0; }
  if (r > rc) { return 0.0; }
  return (-coeffs[2] / powint(r, 2) + coeffs[1] + 3.0 * coeffs[3] * powint(r, 2)) * invnorm;
}

double LddIndicatorSphere::wp2(double r)
{
  if (r < r0) { return 0.0; }
  if (r > rc) { return 0.0; }
  return (2.0 * coeffs[2] / powint(r, 3) + 6.0 * coeffs[3] * r) * invnorm;
}

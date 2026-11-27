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
#include "ldd_indicator_smooth.h"

#include "math_const.h"
#include "math_special.h"
#include "memory.h"

using namespace LAMMPS_NS;
using MathConst::MY_PI;
using MathSpecial::powint;

LddIndicatorSmooth::LddIndicatorSmooth(class LAMMPS *lmp) : LddIndicator(lmp)
{
  n_coeffs = 7;
}

LddIndicatorSmooth::~LddIndicatorSmooth()
{
  if (allocated == 1) { memory->destroy(coeffs); }
  allocated = 0;
}

void LddIndicatorSmooth::allocate()
{
  memory->create(coeffs, n_coeffs, "ldd_indicator:coeffs");
  allocated = 1;
}

void LddIndicatorSmooth::init_coeffs(double a, double b, int dim)
{
  if (!allocated) allocate();
  r0 = a;
  rc = b;
  coeffs[6] = (powint(r0, 5) - powint(rc, 5)) / 120.0 -
      (powint(r0, 3) - powint(rc, 3)) * r0 * rc / 24.0 +
      (r0 - rc) * powint(r0, 2) * powint(rc, 2) / 12.0;
  coeffs[0] =
      (-powint(rc, 5) / 120.0 + r0 * powint(rc, 4) / 24.0 - powint(r0, 2) * powint(rc, 3) / 12.0) /
      coeffs[6];
  coeffs[1] = (powint(r0, 2) * powint(rc, 2) / 4.0) / coeffs[6];
  coeffs[2] = (-rc * r0 * (r0 + rc) / 4.0) / coeffs[6];
  coeffs[3] = ((powint(r0, 2) + 4.0 * r0 * rc + powint(rc, 2)) / 12.0) / coeffs[6];
  coeffs[4] = (-(r0 + rc) / 8.0) / coeffs[6];
  coeffs[5] = 0.05 / coeffs[6];
  switch (dim) {
    case 2:
      norm = 2.0 * MY_PI *
          (powint(r0, 2) / 2.0 + coeffs[5] / 7.0 * (powint(rc, 7) - powint(r0, 7)) +
           coeffs[4] / 6.0 * (powint(rc, 6) - powint(r0, 6)) +
           coeffs[3] / 5.0 * (powint(rc, 5) - powint(r0, 5)) +
           coeffs[2] / 4.0 * (powint(rc, 4) - powint(r0, 4)) +
           coeffs[1] / 3.0 * (powint(rc, 3) - powint(r0, 3)) +
           coeffs[0] / 2.0 * (powint(rc, 2) - powint(r0, 2)));
      break;
    case 3:
      norm = 4.0 * MY_PI *
          (powint(r0, 3) / 3.0 + coeffs[5] / 8.0 * (powint(rc, 8) - powint(r0, 8)) +
           coeffs[4] / 7.0 * (powint(rc, 7) - powint(r0, 7)) +
           coeffs[3] / 6.0 * (powint(rc, 6) - powint(r0, 6)) +
           coeffs[2] / 5.0 * (powint(rc, 5) - powint(r0, 5)) +
           coeffs[1] / 4.0 * (powint(rc, 4) - powint(r0, 4)) +
           coeffs[0] / 3.0 * (powint(rc, 3) - powint(r0, 3)));
      break;
  }
  invnorm = 1.0 / norm;
}

double LddIndicatorSmooth::w(double r)
{
  if (r < r0) { return invnorm; }
  if (r > rc) { return 0.0; }
  return ((coeffs[0] + coeffs[1] * r + coeffs[2] * powint(r, 2) + coeffs[3] * powint(r, 3) +
           coeffs[4] * powint(r, 4) + coeffs[5] * powint(r, 5)) /* coeffs[6]*/) *
      invnorm;
}

double LddIndicatorSmooth::wp(double r)
{
  if (r < r0) { return 0.0; }
  if (r > rc) { return 0.0; }
  return ((coeffs[1] + 2.0 * coeffs[2] * r + 3.0 * coeffs[3] * powint(r, 2) +
           4.0 * coeffs[4] * powint(r, 3) + 5.0 * coeffs[5] * powint(r, 4)) /* coeffs[6]*/) *
      invnorm;
}

double LddIndicatorSmooth::wp2(double r)
{
  if (r < r0) { return 0.0; }
  if (r > rc) { return 0.0; }
  return ((2.0 * coeffs[2] + 6.0 * coeffs[3] * r + 12.0 * coeffs[4] * powint(r, 2) +
           20.0 * coeffs[5] * powint(r, 3)) /* coeffs[6]*/) *
      invnorm;
}

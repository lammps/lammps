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
#include "ldd_indicator_shell.h"

#include "math_const.h"
#include "math_special.h"
#include "memory.h"

using namespace LAMMPS_NS;
using MathConst::MY_PI;
using MathSpecial::powint;

LddIndicatorShell::LddIndicatorShell(class LAMMPS *lmp) : LddIndicator(lmp)
{
  n_coeffs = 7;
}

LddIndicatorShell::~LddIndicatorShell()
{
  if (allocated == 1) memory->destroy(coeffs);
  allocated = 0;
}

void LddIndicatorShell::allocate()
{
  memory->create(coeffs, n_coeffs, "ldd_indicator:coeffs");
  allocated = 1;
}

void LddIndicatorShell::init_coeffs(double a, double b, int dim)
{
  if (!allocated) allocate();
  r0 = a;
  rc = b;
  double denom = powint((1.0 - powint(r0 / rc, 2)), 3);
  coeffs[0] = (1.0 - 3.0 * powint(r0 / rc, 2)) / denom;
  coeffs[2] = 6.0 * powint(r0, 2) / powint(rc, 4) / denom;
  coeffs[4] = -3.0 * (1.0 + powint(r0 / rc, 2)) / powint(rc, 4) / denom;
  coeffs[6] = 2.0 / powint(rc, 6) / denom;
  switch (dim) {
    case 2:
      norm = 2.0 * MY_PI *
          (powint(r0, 2) / 2 + coeffs[0] / 2.0 * (powint(rc, 2) - powint(r0, 2)) +
           coeffs[2] / 4.0 * (powint(rc, 4) - powint(r0, 4)) +
           coeffs[4] / 6.0 * (powint(rc, 6) - powint(r0, 6)) +
           coeffs[6] / 8.0 * (powint(rc, 8) - powint(r0, 8)));
      break;
    case 3:
      norm = 4.0 * MY_PI *
          (powint(r0, 3) / 3 + coeffs[0] / 3.0 * (powint(rc, 3) - powint(r0, 3)) +
           coeffs[2] / 5.0 * (powint(rc, 5) - powint(r0, 5)) +
           coeffs[4] / 7.0 * (powint(rc, 7) - powint(r0, 7)) +
           coeffs[6] / 9.0 * (powint(rc, 9) - powint(r0, 9)));
      break;
  }
  invnorm = 1.0 / norm;
}

double LddIndicatorShell::w(double r)
{
  if (r < r0) { return invnorm; }
  if (r > rc) { return 0.0; }
  return (coeffs[0] + coeffs[2] * powint(r, 2) + coeffs[4] * powint(r, 4) +
          coeffs[6] * powint(r, 6)) *
      invnorm;
}

double LddIndicatorShell::wp(double r)
{
  if (r < r0) { return 0.0; }
  if (r > rc) { return 0.0; }
  return (2.0 * coeffs[2] * r + 4.0 * coeffs[4] * powint(r, 3) + 6.0 * coeffs[6] * powint(r, 5)) *
      invnorm;
}

double LddIndicatorShell::wp2(double r)
{
  if (r < r0) { return 0.0; }
  if (r > rc) { return 0.0; }
  return (2.0 * coeffs[2] + 12.0 * coeffs[4] * powint(r, 2) + 30.0 * coeffs[6] * powint(r, 4)) *
      invnorm;
}

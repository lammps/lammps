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
#include "ldd_indicator_dpd.h"

#include "math_const.h"
#include "math_special.h"
#include "memory.h"

using namespace LAMMPS_NS;
using MathConst::MY_PI;
using MathSpecial::powint;

LddIndicatorDpd::LddIndicatorDpd(class LAMMPS *lmp) : LddIndicator(lmp)
{
  n_coeffs = 3;
}

LddIndicatorDpd::~LddIndicatorDpd()
{
  if (allocated == 1) memory->destroy(coeffs);

  allocated = 0;
}

void LddIndicatorDpd::allocate()
{
  memory->create(coeffs, n_coeffs, "ldd_indicator:coeffs");
  allocated = 1;
}

void LddIndicatorDpd::init_coeffs(double a, double b, int dim)
{
  if (!allocated) allocate();
  r0 = a;
  rc = b;
  coeffs[0] = 1.0;
  coeffs[1] = -2.0 / rc;
  coeffs[2] = 1.0 / powint(rc, 2);
  switch (dim) {
    case 2:
      norm = MY_PI * rc * rc / 6.0;
      break;
    case 3:
      norm = 2.0 * MY_PI * powint(rc, 3) / 15.0;
      break;
  }
  invnorm = 1.0 / norm;
}

double LddIndicatorDpd::w(double r)
{
  if (r > rc) { return 0.0; }
  return (coeffs[0] + coeffs[1] * r + coeffs[2] * powint(r, 2)) * invnorm;
}

double LddIndicatorDpd::wp(double r)
{
  if (r > rc) { return 0.0; }
  return (coeffs[1] + 2.0 * coeffs[2] * r) * invnorm;
}

double LddIndicatorDpd::wp2(double r)
{
  if (r > rc) { return 0.0; }
  return (2.0 * coeffs[2]) * invnorm;
}

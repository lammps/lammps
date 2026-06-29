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
#include "ldd_indicator.h"

using namespace LAMMPS_NS;

LddIndicator::LddIndicator(class LAMMPS * lmp) : Pointers(lmp)
{
  r0 = rc = rs = rb = norm = invnorm = 0.0;
  n_coeffs = allocated = 0;
  coeffs = nullptr;
}

LddIndicator::~LddIndicator()
{
  n_coeffs = 0;
  r0 = rc = rs = rb = norm = invnorm = 0.0;
}

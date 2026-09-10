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

#include "pair_zero_coul.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairZeroCoul::PairZeroCoul(LAMMPS *lmp) : PairZero(lmp)
{
  // claim compatibility with kspace styles.  dispersionflag and tip4pflag
  // must remain unset: those are checked in both directions by
  // KSpace::pair_check(), so setting them would make this pair style
  // incompatible with plain Coulomb kspace styles.

  ewaldflag = pppmflag = msmflag = dipoleflag = spinflag = 1;

  // unlike the parent style, do not offer r-RESPA inner/middle/outer
  // support: the Coulomb pair styles this style stands in for do not
  // support it either

  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void *PairZeroCoul::extract(const char *str, int &dim)
{
  // kspace styles obtain their real-space Coulomb cutoff through this

  if (strcmp(str, "cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_global;
  }
  return nullptr;
}

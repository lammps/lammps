// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: William McDoniel (RWTH Aachen University)
------------------------------------------------------------------------- */

#include "pair_lj_long_coul_long_intel.h"

#include "fix_intel.h"
#include "modify.h"
#include "neigh_request.h"
#include "suffix.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJLongCoulLongIntel::PairLJLongCoulLongIntel(LAMMPS *lmp) :
  PairLJLongCoulLong(lmp)
{
  suffix_flag |= Suffix::INTEL;
  respa_enable = 0;
  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairLJLongCoulLongIntel::init_style()
{
  PairLJLongCoulLong::init_style();

  auto fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");

  fix->pair_init_check();
}

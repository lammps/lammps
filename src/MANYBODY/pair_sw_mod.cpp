// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing authors: Jin-Wu Jiang (Shanghai U) and Wengen Ouyang (Wuhan U)
------------------------------------------------------------------------- */

#include "pair_sw_mod.h"

#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSWMOD::PairSWMOD(LAMMPS *lmp) : PairSW(lmp)
{
  sw_mod_flag = 1;
  delta1 = 0.25;
  delta2 = 0.35;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSWMOD::settings(int narg, char **arg)
{
  // process optional keywords (and not (yet) optional keywords from parent class).
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"maxdelcs") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "pair_style sw/mod", error);
      delta1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      delta2 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;
      if ((delta1 < 0.0) || (delta1 > 1.0) || (delta2 < 0.0) || (delta2 > 1.0) || (delta1 > delta2))
        error->all(FLERR, "Out of range value(s) for pair style sw/mod maxdelcs keyword");
    } else error->all(FLERR, "Illegal pair_style sw/mod keyword: {}", arg[iarg]);
  }
}

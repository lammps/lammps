/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (ICTP, Italy)
------------------------------------------------------------------------- */

#include <cstring>
#include "fix_cac_flux_check.h"
#include "atom.h"
#include "error.h"
#include "output.h"
#include "update.h"
#include "dump.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

CACFixFluxCheck::CACFixFluxCheck(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

}

/* ---------------------------------------------------------------------- */

CACFixFluxCheck::~CACFixFluxCheck()
{

}

/* ---------------------------------------------------------------------- */

int CACFixFluxCheck::setmask()
{
  return PRE_FORCE;
}

/* ---------------------------------------------------------------------- */

void CACFixFluxCheck::init()
{
  //check if CAC atom style is active
  if (!atom->CAC_flag) error->all(FLERR,"CAC fix styles require a CAC atom style");
}

/* ---------------------------------------------------------------------- */

void CACFixFluxCheck::setup_pre_force(int x)
{
  atom->flux_compute = flux_compute = 1;;
}

/* ---------------------------------------------------------------------- */

void CACFixFluxCheck::pre_force(int x)
{
  int ntimestep = update->ntimestep;
  dump = output->dump;
  flux_compute = 0;
  next_dump_any = output->next_dump_any;
  next_dump = output->next_dump;
  last_dump = output->last_dump;
  ndump = output->ndump;

  if (next_dump_any == ntimestep) {
    for (int idump = 0; idump < ndump; idump++) {
      if (next_dump[idump] == ntimestep) {
        if (last_dump[idump] != ntimestep) {
          if (strcmp(dump[idump]->style, "cac/flux") == 0 || strcmp(dump[idump]->style, "cac/atom/flux") == 0)
            flux_compute = 1;
        }
      }
    }
  }

  atom->flux_compute = flux_compute;
}

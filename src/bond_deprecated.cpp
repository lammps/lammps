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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "bond_deprecated.h"

#include "bond_hybrid.h"
#include "comm.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void BondDeprecated::settings(int, char **)
{
  std::string my_style = force->bond_style;

  // hybrid substyles are created in BondHybrid::settings(), so when this is
  // called, our style was just added at the end of the list of substyles

  if (utils::strmatch(my_style,"^hybrid")) {
    BondHybrid *hybrid = (BondHybrid *)force->bond;
    my_style = hybrid->keywords[hybrid->nstyles];
  }

  if (my_style == "DEPRECATED") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp,"\nBond style 'DEPRECATED' is a dummy style\n\n");
    return;
  }
  error->all(FLERR,"This bond style is no longer available");
}



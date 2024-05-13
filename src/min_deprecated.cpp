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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "min_deprecated.h"

#include "comm.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MinDeprecated::MinDeprecated(LAMMPS *_lmp) : Min(_lmp)
{
  std::string my_style = update->minimize_style;

  if (my_style == "DEPRECATED") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp, "\nMinimize style 'DEPRECATED' is a dummy style\n\n");
    return;
  }

  if (my_style == "fire/old")
    error->all(FLERR,
               "Minimize style 'fire/old' has been removed from LAMMPS after "
               "the 22 December 2022 version.\nERROR: Please use 'min_style fire'");

  error->all(FLERR, "This minimize style is no longer available");
}

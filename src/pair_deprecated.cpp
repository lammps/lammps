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

#include "pair_deprecated.h"

#include "comm.h"
#include "error.h"
#include "force.h"
#include "pair_hybrid.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void PairDeprecated::settings(int, char **)
{
  std::string my_style = force->pair_style;

  // hybrid substyles are created in PairHybrid::settings(), so when this is
  // called, our style was just added at the end of the list of substyles

  if (utils::strmatch(my_style, "^hybrid")) {
    auto hybrid = dynamic_cast<PairHybrid *>(force->pair);
    my_style = hybrid->keywords[hybrid->nstyles];
  }

  if (my_style == "DEPRECATED") {
    if (lmp->comm->me == 0) utils::logmesg(lmp, "\nPair style 'DEPRECATED' is a dummy style\n\n");
    return;
  }

  if (my_style == "reax") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp,
                     "\nPair style 'reax' has been removed from LAMMPS "
                     "after the 12 December 2018 version\n\n");
  } else if (my_style == "mesont/tpm") {
    if (lmp->comm->me == 0)
      utils::logmesg(lmp,
                     "\nPair style 'mesont/tpm' has been removed from LAMMPS. "
                     "Please use pair style 'mesocnt' instead\n\n");
  }
  error->all(FLERR, "This pair style is no longer available");
}

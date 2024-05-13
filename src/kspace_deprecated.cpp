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

#include "kspace_deprecated.h"

#include "comm.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void KSpaceDeprecated::settings(int, char **)
{
  std::string my_style = force->kspace_style;

  if (my_style == "DEPRECATED") {
    if (lmp->comm->me == 0) utils::logmesg(lmp, "\nKSpace style 'DEPRECATED' is a dummy style\n\n");
    return;
  }
  error->all(FLERR, "This kspace style is no longer available");
}

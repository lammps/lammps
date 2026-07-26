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

#include "fix_pimd_nvt_bosonic.h"

#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixPIMDBNVT::FixPIMDBNVT(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  error->all(FLERR, "Fix pimd/nvt/bosonic is deprecated; use fix pimd/langevin/bosonic instead");
}

/* ---------------------------------------------------------------------- */

int FixPIMDBNVT::setmask()
{
  return 0;
}

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
   Contributing author:  Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "fix_colvars_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

FixColvarsKokkos::FixColvarsKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixColvars(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

void FixColvarsKokkos::post_force(int vflag)
{
  atomKK->sync(Host,X_MASK|F_MASK|TAG_MASK|IMAGE_MASK);
  FixColvars::post_force(vflag);
  atomKK->modified(Host,F_MASK);
}

/* ---------------------------------------------------------------------- */
void FixColvarsKokkos::end_of_step()
{
  if (store_forces) {
    atomKK->sync(Host,F_MASK|TAG_MASK);
    FixColvars::end_of_step();
  }
}

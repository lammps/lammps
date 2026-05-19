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

#include "fix_deform_pressure_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixDeformPressureKokkos::FixDeformPressureKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDeformPressure(lmp, narg, arg)
{
  kokkosable = 1;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ----------------------------------------------------------------------
   box flipped on previous step: sync all per-atom data to host first
------------------------------------------------------------------------- */

void FixDeformPressureKokkos::pre_exchange()
{
  atomKK->sync(Host, ALL_MASK);
  FixDeformPressure::pre_exchange();
  atomKK->modified(Host, ALL_MASK);
}

/* ---------------------------------------------------------------------- */

void FixDeformPressureKokkos::update_box()
{
  // update_box() may remap atom positions when remapflag == X_REMAP
  if (remapflag == Domain::X_REMAP && rfix.size() > 0)
    atomKK->sync(Host, ALL_MASK);

  FixDeformPressure::update_box();

  if (remapflag == Domain::X_REMAP && rfix.size() > 0)
    atomKK->modified(Host, ALL_MASK);
}

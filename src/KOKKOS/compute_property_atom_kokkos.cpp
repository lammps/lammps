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

#include "compute_property_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

ComputePropertyAtomKokkos:: ComputePropertyAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputePropertyAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyAtomKokkos::compute_peratom()
{
  atomKK->sync(Host, ALL_MASK);
  ComputePropertyAtom::compute_peratom();
}

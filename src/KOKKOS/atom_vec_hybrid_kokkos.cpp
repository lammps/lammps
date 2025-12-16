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

#include "atom_vec_hybrid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecHybridKokkos::AtomVecHybridKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecHybrid(lmp)
{
}

/* ----------------------------------------------------------------------
   process field strings to initialize data structs for all other methods
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::init()
{
  AtomVecHybrid::init();

  set_atom_masks();
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::grow(int n)
{
  for (int k = 0; k < nstyles; k++) styles[k]->grow(n);
  nmax = atomKK->k_x.view_host().extent(0);

  tag = atom->tag;
  type = atom->type;
  mask = atom->mask;
  image = atom->image;
  x = atom->x;
  v = atom->v;
  f = atom->f;
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  for (int k = 0; k < nstyles; k++)
    (dynamic_cast<AtomVecKokkos*>(styles[k]))->sort_kokkos(Sorter);
}

// TODO: move dynamic_cast into init

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::sync(ExecutionSpace space, uint64_t h_mask)
{
  for (int k = 0; k < nstyles; k++) (dynamic_cast<AtomVecKokkos*>(styles[k]))->sync(space,h_mask);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::sync_pinned(ExecutionSpace space, uint64_t h_mask, int async_flag)
{
  for (int k = 0; k < nstyles; k++) (dynamic_cast<AtomVecKokkos*>(styles[k]))->sync_pinned(space,h_mask,async_flag);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::modified(ExecutionSpace space, uint64_t h_mask)
{
  for (int k = 0; k < nstyles; k++) (dynamic_cast<AtomVecKokkos*>(styles[k]))->modified(space,h_mask);
}

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

#include "fix_oxdna_prime_neighs_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neighbor_kokkos.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixOxdnaPrimeNeighsKokkos<DeviceType>::FixOxdnaPrimeNeighsKokkos(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = TAG_MASK | CG_DNA_MASK;
  datamask_modify = EMPTY_MASK;

  nbondlist = 0;
  last_atom_ncalls = -1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixOxdnaPrimeNeighsKokkos<DeviceType>::~FixOxdnaPrimeNeighsKokkos() = default;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::init()
{
  // No neighbor list requested: all lookups are computed on demand by the
  // styles that use them.

  // atoms may have been reordered since the last run
  last_atom_ncalls = -1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixOxdnaPrimeNeighsKokkos<DeviceType>::setmask()
{
  // all lookups are computed on demand by the styles that use them

  return 0;
}

/* ----------------------------------------------------------------------
   Called by the bond (fene) and pair (stk) styles to precompute the
   prime-neighbor lookups for the bond list that is current at the time
   of the call.  Each caller passes its own output View, since under
   bond style hybrid the bond style is handed a sub-style bond list
   while the stk pair styles loop over the full bond list.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::compute_prime_neighs_bond(
  typename AT::t_int_1d_4 &d_prime_neighs)
{
  neighborKK->k_bondlist.template sync<DeviceType>();
  bondlist = neighborKK->k_bondlist.view<DeviceType>();
  nbondlist = neighborKK->nbondlist;

  if (nbondlist > d_prime_neighs.extent_int(0))
    MemKK::realloc_kokkos(d_prime_neighs, "prime_neighs:prime_neighs_bond", nbondlist);
  d_prime_neighs_bond = d_prime_neighs;

  atomKK->sync(execution_space, datamask_read);
  tag = atomKK->k_tag.view<DeviceType>();
  id5p = atomKK->k_id5p.view<DeviceType>();
  id3p = atomKK->k_id3p.view<DeviceType>();

  map_style = atom->map_style;
  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }

  copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType, TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsBond>(0, nbondlist),
    *this);
  copymode = 0;
}

/* ----------------------------------------------------------------------
   Called by the excv and oxdna3/xstk pair styles from their compute()
   to precompute the local indices of the 3' and 5' neighbors of every
   local and ghost atom.  Shared by both styles, so only rebuilt once per
   reneighbor.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::compute_prime_neighs_atom()
{
  if (last_atom_ncalls == neighbor->ncalls) return;
  last_atom_ncalls = neighbor->ncalls;

  const int nall = atom->nlocal + atom->nghost;

  if (nall > d_prime_neighs_atom.extent_int(0))
    MemKK::realloc_kokkos(d_prime_neighs_atom, "prime_neighs:prime_neighs_atom", nall, 2);

  atomKK->sync(execution_space, datamask_read);
  tag = atomKK->k_tag.view<DeviceType>();
  id5p = atomKK->k_id5p.view<DeviceType>();
  id3p = atomKK->k_id3p.view<DeviceType>();

  map_style = atom->map_style;
  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }

  copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType, TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsAtom>(0, nall),
    *this);
  copymode = 0;
}

/* ----------------------------------------------------------------------
   Loop through the bondlist and precompute the atom mapping for
   the 3' and 5' neighbors of each bonded pair. This is the KOKKOS
   equivalent of "atom->map(id{3/5}p[{a/b}])" in the CPU code.
   These indexes are then used directly within the main compute loop.
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixOxdnaPrimeNeighsKokkos<DeviceType>::operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsBond,
                                                       const int &in) const
{
  // Bondlist contains local atom indices (can be >= nlocal for ghosts).
  // [k/d]_bondlist already has KOKKOS 'closest_image' applied, so we can use these directly.
  int a = bondlist(in,0);
  int b = bondlist(in,1);

  // Directionality test: a -> b must be 3' -> 5'
  int atom_a = a;
  int atom_b = b;
  if (tag(b) != id5p(a)) {
    atom_a = b;
    atom_b = a;
  }

  d_prime_neighs_bond(in,0) = atom_a;
  d_prime_neighs_bond(in,1) = atom_b;

  // Look up local indices of the 3'/5' tetramer-context neighbors.
  // These are only used for type() lookup in the main compute loop,
  // so map_kokkos (tag -> local index) is sufficient; no closest_image needed.
  //
  // We break the oxDNA: 3'neighbor(a) - a - b - 5'neighbor(b) convention here.
  // Instead, we have: a, b, 3'neighbor(a), 5'neighbor(b) - this is the order that
  // they are actually accessed in the main compute loop.
  //
  int id3p_local = -1;
  const tagint id3p_tag = id3p(atom_a);
  int mapped = -1;
  if (id3p_tag != -1) {
    if (map_style == Atom::MAP_ARRAY) {
      const auto map_array = k_map_array.view<DeviceType>();
      if (id3p_tag >= 0 && id3p_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id3p_tag);
    } else if (map_style == Atom::MAP_HASH) {
      mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id3p_tag, k_map_hash);
    }
    if (mapped >= 0) id3p_local = mapped;
    // If a present id3p tag ever fails to map locally, the downstream stk/fene
    // code will treat it like a terminal neighbor (tetramer type 0). In valid
    // oxDNA runs this path is expected to be unreachable because bonded context
    // atoms needed for tetramer typing should also be present as local or ghost.
  }
  d_prime_neighs_bond(in,2) = id3p_local;

  int id5p_local = -1;
  const tagint id5p_tag = id5p(atom_b);
  if (id5p_tag != -1) {
    mapped = -1;
    if (map_style == Atom::MAP_ARRAY) {
      const auto map_array = k_map_array.view<DeviceType>();
      if (id5p_tag >= 0 && id5p_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id5p_tag);
    } else if (map_style == Atom::MAP_HASH) {
      mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id5p_tag, k_map_hash);
    }
    if (mapped >= 0) id5p_local = mapped;
    // Same assumption as above for id5p: a map miss falls back to the terminal
    // tetramer type sentinel but is expected to be unreachable in normal runs.
  }
  d_prime_neighs_bond(in,3) = id5p_local;
}

/* ----------------------------------------------------------------------
   Precompute the local indices of the 3' and 5' neighbors of atom i.
   This is the KOKKOS equivalent of "atom->map(id{3/5}p[i])" in the CPU
   code; -1 marks a missing (terminal) neighbor.
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixOxdnaPrimeNeighsKokkos<DeviceType>::operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsAtom,
                                                       const int &i) const
{
  int mapped = -1;
  const tagint id3p_tag = id3p(i);
  if (id3p_tag != -1) {
    if (map_style == Atom::MAP_ARRAY) {
      const auto map_array = k_map_array.view<DeviceType>();
      if (id3p_tag >= 0 && id3p_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id3p_tag);
    } else if (map_style == Atom::MAP_HASH) {
      mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id3p_tag, k_map_hash);
    }
  }
  d_prime_neighs_atom(i,0) = mapped;

  mapped = -1;
  const tagint id5p_tag = id5p(i);
  if (id5p_tag != -1) {
    if (map_style == Atom::MAP_ARRAY) {
      const auto map_array = k_map_array.view<DeviceType>();
      if (id5p_tag >= 0 && id5p_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id5p_tag);
    } else if (map_style == Atom::MAP_HASH) {
      mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id5p_tag, k_map_hash);
    }
  }
  d_prime_neighs_atom(i,1) = mapped;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixOxdnaPrimeNeighsKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixOxdnaPrimeNeighsKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS

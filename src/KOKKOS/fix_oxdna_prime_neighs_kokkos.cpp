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
#include "fix_oxdna_npair_kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list_kokkos.h"
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
  anum = 0;
  npairlist = 0;
  fix_oxdna_npairKK = nullptr;
  last_precompute_lastcall = -1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixOxdnaPrimeNeighsKokkos<DeviceType>::~FixOxdnaPrimeNeighsKokkos() = default;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::init()
{
  // No neighbor list requested: the pair style supplies its own list to
  // compute_prime_neighs_pair().  Bond precompute uses the global bondlist.
  // NOTE: We do not hard-require OXDNA/NPAIR/kk at init time.  Some style
  // combinations initialize PRIME_NEIGHS before NPAIR.  We resolve NPAIR
  // lazily in compute_prime_neighs_oxdna3_xstk(), which is the only path
  // that needs it.
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixOxdnaPrimeNeighsKokkos<DeviceType>::setmask()
{
  int mask = 0;
  mask |= MIN_PRE_FORCE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::min_setup_pre_force(int vflag)
{
  min_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::min_pre_force(int /*vflag*/)
{
  if (neighbor->lastcall != last_precompute_lastcall) {
    compute_prime_neighs_bond();
    last_precompute_lastcall = neighbor->lastcall;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::pre_force(int /*vflag*/)
{
  if (neighbor->lastcall != last_precompute_lastcall) {
    compute_prime_neighs_bond();
    last_precompute_lastcall = neighbor->lastcall;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::compute_prime_neighs_bond()
{
  // Bond precompute only. Pair precompute is driven by the pair style via
  // compute_prime_neighs_pair() so that it always uses the pair's own list.
  neighborKK->k_bondlist.template sync<DeviceType>();
  bondlist = neighborKK->k_bondlist.view<DeviceType>();
  nbondlist = neighborKK->nbondlist;

  if (nbondlist > d_prime_neighs_bond.extent_int(0)) {
    MemKK::realloc_kokkos(d_prime_neighs_bond, "prime_neighs:prime_neighs_bond", nbondlist);
  }

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
   Called by the pair style from its compute() to precompute prime-neighbor
   map lookups for the pair neighbor list.  Uses the pair's own list so
   that row index "ib" in d_prime_neighs_pair(a,ib,...) always matches the
   ib-th entry in the pair's d_neighbors(a,ib).
------------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::compute_prime_neighs_pair(NeighList *neigh_list)
{
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(neigh_list);
  d_neighbors = k_list->d_neighbors;
  anum = neigh_list->inum;
  d_alist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  const int maxneigh = d_neighbors.extent(1);

  if (anum > d_prime_neighs_pair.extent_int(0) || maxneigh > d_prime_neighs_pair.extent_int(1)) {
    MemKK::realloc_kokkos(k_prime_neighs_pair, "prime_neighs:prime_neighs_pair", anum, maxneigh, 4);
    d_prime_neighs_pair = k_prime_neighs_pair.template view<DeviceType>();
  }

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
    Kokkos::RangePolicy<DeviceType, TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsPair>(0, anum),
    *this);
  copymode = 0;
}

/* ----------------------------------------------------------------------
   Called (currently) only from oxdna3/xstk from its compute() to
   precompute prime-neighbor map lookups for its npair neighbor list,
   which comes from fix_oxdna_npair_kokkos*.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaPrimeNeighsKokkos<DeviceType>::compute_prime_neighs_oxdna3_xstk(NeighList *neigh_list)
{
  (void) neigh_list;

  if (!fix_oxdna_npairKK) {
    auto npair_fixes = modify->get_fix_by_style("^OXDNA/NPAIR/kk");
    for (auto *fixptr : npair_fixes) {
      auto *typed = dynamic_cast<FixOxdnaNpairKokkos<DeviceType> *>(fixptr);
      if (typed) {
        fix_oxdna_npairKK = typed;
        break;
      }
    }
  }

  if (!fix_oxdna_npairKK) {
    error->all(FLERR, "FixOxdnaPrimeNeighsKokkos::compute_prime_neighs_oxdna3_xstk() "
               "called but no matching OXDNA/NPAIR/kk fix is available");
  }

  npairlist = fix_oxdna_npairKK->screened_pair_count;
  pairlist = fix_oxdna_npairKK->k_pairs_screened.template view<DeviceType>();

  if (npairlist > d_prime_neighs_oxdna3_xstk.extent_int(0)) {
    MemKK::realloc_kokkos(k_prime_neighs_oxdna3_xstk,
                          "prime_neighs:prime_neighs_oxdna3_xstk", npairlist, 4);
    d_prime_neighs_oxdna3_xstk = k_prime_neighs_oxdna3_xstk.template view<DeviceType>();
  }

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
    Kokkos::RangePolicy<DeviceType, TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsOxdna3Xstk>(0, npairlist),
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
   Loop through the neighbor list and precompute the atom mapping for
   the 3' and 5' neighbors of each pair. This is the KOKKOS
   equivalent of "atom->map(id{3/5}p[{a/b}])" in the CPU code.
   These indexes are then used directly within the main compute loop.
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixOxdnaPrimeNeighsKokkos<DeviceType>::operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsPair,
                                                       const int &in) const
{
  // excv runs thread-per-particle loop. Find one same-strand nearest-neighbor
  // partner b (if present) and precompute the map() lookups used in the vanilla
  // excv base-base topology branch.
  const int a = d_alist(in);
  const int jnum = d_numneigh(a);

  for (int jj = 0; jj < jnum; jj++) {
    const int btry = d_neighbors(a,jj) & NEIGHMASK;

    int mapped = -1; // default to -1 for missing neighbor, which is treated as a terminal neighbor in downstream compute path.

    const tagint id3p_a_tag = id3p(a);
    if (id3p_a_tag != -1) {
      if (map_style == Atom::MAP_ARRAY) {
        const auto map_array = k_map_array.view<DeviceType>();
        if (id3p_a_tag >= 0 && id3p_a_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id3p_a_tag);
      } else if (map_style == Atom::MAP_HASH) {
        mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id3p_a_tag, k_map_hash);
      }
    }
    d_prime_neighs_pair(a,jj,0) = mapped;

    mapped = -1;
    const tagint id5p_b_tag = id5p(btry);
    if (id5p_b_tag != -1) {
      if (map_style == Atom::MAP_ARRAY) {
        const auto map_array = k_map_array.view<DeviceType>();
        if (id5p_b_tag >= 0 && id5p_b_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id5p_b_tag);
      } else if (map_style == Atom::MAP_HASH) {
        mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id5p_b_tag, k_map_hash);
      }
    }
    d_prime_neighs_pair(a,jj,1) = mapped;

    mapped = -1;
    const tagint id3p_b_tag = id3p(btry);
    if (id3p_b_tag != -1) {
      if (map_style == Atom::MAP_ARRAY) {
        const auto map_array = k_map_array.view<DeviceType>();
        if (id3p_b_tag >= 0 && id3p_b_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id3p_b_tag);
      } else if (map_style == Atom::MAP_HASH) {
        mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id3p_b_tag, k_map_hash);
      }
    }
    d_prime_neighs_pair(a,jj,2) = mapped;

    mapped = -1;
    const tagint id5p_a_tag = id5p(a);
    if (id5p_a_tag != -1) {
      if (map_style == Atom::MAP_ARRAY) {
        const auto map_array = k_map_array.view<DeviceType>();
        if (id5p_a_tag >= 0 && id5p_a_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id5p_a_tag);
      } else if (map_style == Atom::MAP_HASH) {
        mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id5p_a_tag, k_map_hash);
      }
    }
    d_prime_neighs_pair(a,jj,3) = mapped;
  }
}

/* ----------------------------------------------------------------------
   Loop through custom oxdna npair neighbor list and precompute the atom mapping for
   the 3' and 5' neighbors of each pair. This is the KOKKOS equivalent of
   "atom->map(id{3/5}p[{a/b}]) in the CPU code. These indexes are then used
   directly within the main compute loop.
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixOxdnaPrimeNeighsKokkos<DeviceType>::operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsOxdna3Xstk,
                                                       const int &ipair) const
{
  // Direct packed pair lookup: high 32 bits = a, low 32 bits = b.
  const uint64_t pair = pairlist(ipair);
  // "pair >> 32" shifts the pair to the right by 32 bits, so the upper 32 bits
  // becomes the lower 32 bits to recover the atom-a index.
  const int a = static_cast<int>(pair >> 32);
  // "pair & 0xffffffffu" keeps only the lower 32 bits to recover the atom-b index.
  int b = static_cast<int>(pair & 0xffffffffu);
  b &= NEIGHMASK;

  int mapped = -1;

  // id3p[a]
  const tagint id3p_a_tag = id3p(a);
  if (id3p_a_tag != -1) {
    if (map_style == Atom::MAP_ARRAY) {
      const auto map_array = k_map_array.view<DeviceType>();
      if (id3p_a_tag >= 0 && id3p_a_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id3p_a_tag);
    } else if (map_style == Atom::MAP_HASH) {
      mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id3p_a_tag, k_map_hash);
    }
  }
  d_prime_neighs_oxdna3_xstk(ipair,0) = mapped;

  // id5p[a]
  mapped = -1;
  const tagint id5p_a_tag = id5p(a);
  if (id5p_a_tag != -1) {
    if (map_style == Atom::MAP_ARRAY) {
      const auto map_array = k_map_array.view<DeviceType>();
      if (id5p_a_tag >= 0 && id5p_a_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id5p_a_tag);
    } else if (map_style == Atom::MAP_HASH) {
      mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id5p_a_tag, k_map_hash);
    }
  }
  d_prime_neighs_oxdna3_xstk(ipair,1) = mapped;

  // id3p[b]
  mapped = -1;
  const tagint id3p_b_tag = id3p(b);
  if (id3p_b_tag != -1) {
    if (map_style == Atom::MAP_ARRAY) {
      const auto map_array = k_map_array.view<DeviceType>();
      if (id3p_b_tag >= 0 && id3p_b_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id3p_b_tag);
    } else if (map_style == Atom::MAP_HASH) {
      mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id3p_b_tag, k_map_hash);
    }
  }
  d_prime_neighs_oxdna3_xstk(ipair,2) = mapped;

  // id5p[b]
  mapped = -1;
  const tagint id5p_b_tag = id5p(b);
  if (id5p_b_tag != -1) {
    if (map_style == Atom::MAP_ARRAY) {
      const auto map_array = k_map_array.view<DeviceType>();
      if (id5p_b_tag >= 0 && id5p_b_tag < static_cast<tagint>(map_array.extent(0))) mapped = map_array(id5p_b_tag);
    } else if (map_style == Atom::MAP_HASH) {
      mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(id5p_b_tag, k_map_hash);
    }
  }
  d_prime_neighs_oxdna3_xstk(ipair,3) = mapped;

}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixOxdnaPrimeNeighsKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixOxdnaPrimeNeighsKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS

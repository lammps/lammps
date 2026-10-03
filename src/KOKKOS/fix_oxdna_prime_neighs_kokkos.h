/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(OXDNA/PRIME_NEIGHS/kk,FixOxdnaPrimeNeighsKokkos<LMPDeviceType>);
FixStyle(OXDNA/PRIME_NEIGHS/kk/device,FixOxdnaPrimeNeighsKokkos<LMPDeviceType>);
FixStyle(OXDNA/PRIME_NEIGHS/kk/host,FixOxdnaPrimeNeighsKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_OXDNA_PRIME_NEIGHS_KOKKOS_H
#define LMP_FIX_OXDNA_PRIME_NEIGHS_KOKKOS_H

#include "fix.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsBond {}; // fene and stk

struct TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsAtom {}; // excv and oxdna3/xstk

template<class DeviceType>
class FixOxdnaPrimeNeighsKokkos : public Fix {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixOxdnaPrimeNeighsKokkos(class LAMMPS *, int, char **);
  ~FixOxdnaPrimeNeighsKokkos() override;

  void init() override;
  int setmask() override;

  // ------ For PrimeNeighBond (fene and stk)
  // 0-3 : atom a, atom b, id3p[a], id5p[b] for each bond of the current bond list.
  // As per their order of being called in fene and stk compute.
  // The caller owns the output View (grown as needed).
  void compute_prime_neighs_bond(typename AT::t_int_1d_4 &d_prime_neighs);
  // ------ For PrimeNeighAtom (excv and oxdna3/xstk)
  // 0-1 : local index of id3p[i], id5p[i] (-1 if none) for each local and ghost atom i.
  // Recomputed at most once per reneighbor; callers must re-fetch the View afterwards.
  typename AT::t_int_2d d_prime_neighs_atom;
  void compute_prime_neighs_atom();

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsBond, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsAtom, const int &) const;

 private:
  class NeighborKokkos *neighborKK;

  typename AT::t_tagint_1d tag;
  typename AT::t_tagint_1d id5p;
  typename AT::t_tagint_1d id3p;
  // For PrimeNeighBond
  int nbondlist;
  typename AT::t_int_2d_lr bondlist;
  typename AT::t_int_1d_4 d_prime_neighs_bond;
  // For PrimeNeighAtom
  bigint last_atom_ncalls;

  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;
};

}    // namespace LAMMPS_NS
#endif
#endif

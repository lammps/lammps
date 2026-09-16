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

template<class DeviceType>
class FixOxdnaNpairKokkos;

struct TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsBond {}; // fene and stk

struct TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsPair {}; // excv

struct TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsOxdna3Xstk {}; // oxdna3/xstk

template<class DeviceType>
class FixOxdnaPrimeNeighsKokkos : public Fix {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixOxdnaPrimeNeighsKokkos(class LAMMPS *, int, char **);
  ~FixOxdnaPrimeNeighsKokkos() override;

  void init() override;
  int setmask() override;
  void min_setup_pre_force(int);
  void min_pre_force(int) override;
  void setup_pre_force(int) override;
  void pre_force(int) override;

  // ------ For PrimeNeighBond (fene and stk)
  // 0-3 : atom a, atom b, id3p[a], id5p[b] for each bond.
  // As per their order of being called in fene and stk compute.
  typename AT::t_int_1d_4 d_prime_neighs_bond;
  void compute_prime_neighs_bond();
  // ------ For PrimeNeighPair (excv)
  // 0-3 : id3p[a], id5p[b], id3p[b], id5p[a] for each pair.
  // As per their order of being called in excv compute.
  // Layout matches the native neighlist walk: d_prime_neighs_pair(a,ib,0-3).
  // Populated by compute_prime_neighs_pair(), called by the pair style from
  // its compute() using the pair's own neighbor list.
  DAT::tdual_int_3d k_prime_neighs_pair;
  typename AT::t_int_3d d_prime_neighs_pair;
  void compute_prime_neighs_pair(class NeighList *neigh_list);
  // ------ For PrimeNeighOxdna3Xstk (oxdna3/xstk/kk)
  // 0-3 : id3p[a], id5p[b], id3p[b], id5p[a] for each pair.
  // As per their order of being called in oxdna3/xstk compute.
  // Layout is per screened pair index from fix_oxdna_npair_kokkos:
  // d_prime_neighs_oxdna3_xstk(ipair,0-3), where ipair maps to the packed
  // (a,braw) pair in npair's d_pairs_screened.
  // Populated by compute_prime_neighs_oxdna3_xstk(), called by the pair style
  // from its compute() using the pair's own neighbor list.
  DAT::tdual_int_2d k_prime_neighs_oxdna3_xstk;
  typename AT::t_int_2d d_prime_neighs_oxdna3_xstk;
  void compute_prime_neighs_oxdna3_xstk(class NeighList *neigh_list);

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsBond, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsPair, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsOxdna3Xstk, const int&) const;

 private:
  class NeighborKokkos *neighborKK;

  typename AT::t_tagint_1d tag;
  typename AT::t_tagint_1d id5p;
  typename AT::t_tagint_1d id3p;
  // For PrimeNeighBond
  int nbondlist;
  typename AT::t_int_2d_lr bondlist;
  // For PrimeNeighPair (set in compute_prime_neighs_pair)
  int anum;
  typename AT::t_neighbors_2d_randomread d_neighbors;
  typename AT::t_int_1d_randomread d_alist;
  typename AT::t_int_1d_randomread d_numneigh;
  // For PrimeNeighOxdna3Xstk (set in compute_prime_neighs_oxdna3_xstk)
  int npairlist;
  typename AT::t_uint64_1d pairlist;

  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  bigint last_precompute_lastcall;

  FixOxdnaNpairKokkos<DeviceType> *fix_oxdna_npairKK;    // ptr to OXDNA/NPAIR/kk fix
};

}    // namespace LAMMPS_NS
#endif
#endif

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


struct TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsOxdna3Xstk {}; // oxdna3/xstk

struct TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsAtom {}; // excv

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
  // ------ For PrimeNeighOxdna3Xstk (oxdna3/xstk/kk)
  // 0-3 : id3p[a], id5p[a], id3p[b], id5p[b] for each pair.
  // As per their order of being called in oxdna3/xstk compute.
  // Layout is per screened pair index from fix_oxdna_npair_kokkos:
  // d_prime_neighs_oxdna3_xstk(ipair,0-3), where ipair maps to the packed
  // (a,braw) pair in npair's d_pairs_screened.
  // Populated by compute_prime_neighs_oxdna3_xstk(), called by the pair style
  // from its compute() whenever the screened list was rebuilt.
  DAT::tdual_int_2d k_prime_neighs_oxdna3_xstk;
  typename AT::t_int_2d d_prime_neighs_oxdna3_xstk;
  void compute_prime_neighs_oxdna3_xstk();
  // ------ For PrimeNeighAtom (excv)
  // 0-3 : local index of id3p[i] and id5p[i] (-1 if none), and their atom
  // types (0 if none), for each owned and ghost atom i
  typename AT::t_int_1d_4 d_prime_neighs_atom;
  void compute_prime_neighs_atom();

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsBond, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsOxdna3Xstk, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaPrimeNeighsPrecomputePrimeNeighsAtom, const int&) const;

 private:
  class NeighborKokkos *neighborKK;

  typename AT::t_tagint_1d tag;
  typename AT::t_tagint_1d id5p;
  typename AT::t_tagint_1d id3p;
  // For PrimeNeighBond
  int nbondlist;
  typename AT::t_int_2d_lr bondlist;
  typename AT::t_int_1d_4 d_prime_neighs_bond;
  // For PrimeNeighOxdna3Xstk (set in compute_prime_neighs_oxdna3_xstk)
  int npairlist;
  typename AT::t_uint64_1d pairlist;

  typename AT::t_int_1d type;

  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;
};

}    // namespace LAMMPS_NS
#endif
#endif

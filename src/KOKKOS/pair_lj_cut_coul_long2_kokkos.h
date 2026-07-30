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

/* ----------------------------------------------------------------------
   lj/cut/coul/long2/kk -- experimental hand-written GPU variant of
   lj/cut/coul/long/kk.  Physics is identical; the difference is the
   force-only compute path, which bypasses the generic pair_kokkos.h
   dispatch (PairComputeFunctor / ScatterView / CoulLongTable templating)
   in favor of one lean kernel: a single fused loop with analytical Ewald
   (no coul table lookup), register-resident i-force accumulation, and
   direct atomic j-scatter.  Energy/virial (thermo) steps delegate to the
   proven base-class kernel, so only the hot force-only path is new.

   Scope: GPU, single pair style, half list + newton on (-> HALFTHREAD),
   ntypes <= MAX_TYPES_STACKPARAMS.  Any other configuration transparently
   falls back to PairLJCutCoulLongKokkos::compute().  See kokkos_neigh.md.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/long2/kk,PairLJCutCoulLong2Kokkos<LMPDeviceType>);
PairStyle(lj/cut/coul/long2/kk/device,PairLJCutCoulLong2Kokkos<LMPDeviceType>);
PairStyle(lj/cut/coul/long2/kk/host,PairLJCutCoulLong2Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_COUL_LONG2_KOKKOS_H
#define LMP_PAIR_LJ_CUT_COUL_LONG2_KOKKOS_H

#include "pair_lj_cut_coul_long_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJCutCoulLong2Kokkos : public PairLJCutCoulLongKokkos<DeviceType> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairLJCutCoulLong2Kokkos(class LAMMPS *);

  void compute(int, int) override;

 protected:
  // ---- i-batched union neighbor list (experimental, LMP_LJCL2_UNION=4|8) ----
  // CI consecutive i atoms of the spatially sorted d_ilist share ONE entry per
  // distinct neighbor j.  Entry layout, one int:
  //   bits  0..21  j
  //   bits 22..29  member mask, bit k set iff (i_k, j) is a flat-list entry
  //   bits 30..31  special-bond code (only ever on a single-member entry)
  // The mask makes the walked pair set identical to the flat half list -- no
  // pair computed twice, none lost -- while the per-ENTRY costs (index load,
  // x/q/type gather, atomic f[j] scatter) are amortized over popcount members.
  // See kokkos_neigh.md section 20.
  typename AT::t_neighbors_2d_lr d_unbr;
  typename AT::t_int_1d d_unum;
  DAT::tdual_int_1d k_ustat;    // [0] overflow, [1] max entries/group, [2] entries,
                                // [3] pairs, [4] specials
  int union_ci = -1;            // -1 = not yet probed, 0 = off, else CI
  int union_maxu = 0;
  int union_ngroup = 0;
  int union_stats = 0;
  bigint union_built_step = -1;

  // ---- fused dual-cutoff inner list (experimental, LMP_LJCL2_DUAL=1) ----
  // Half the flat kernel's time goes to skin-shell pairs that fail the force
  // cutoff: at skin 6.0 it walks 755 neighbors/atom of which 185 are inside
  // the 10 A cutoff.  Every inner_every steps (and whenever any atom has moved
  // more than half the inner skin) the kernel walks the master list as usual
  // and, as a by-product of the distance test it already computes, stores the
  // pairs within cutforce+inner_skin into d_inbr.  The other steps walk that
  // list -- ~1/3 the entries, so the removed work is gone from the warp's
  // issue stream rather than merely masked off.  Unlike neigh/prune (sections
  // 14, 19) there is no separate full-traffic prune pass: the marginal cost of
  // a refresh step is one store per kept pair.  See kokkos_neigh.md section 20.
  typename AT::t_neighbors_2d_lr d_inbr;      // inner list, LayoutRight (jj-fast)
  typename AT::t_int_1d d_innum;
  typename AT::t_kkfloat_1d_3_lr d_xhold;     // positions at the last refresh
  Kokkos::View<KK_FLOAT, DeviceType> d_maxsq; // max squared displacement since then
  DAT::tdual_int_1d k_ictl;   // [0] refresh flag (device-side control, no host
                              //     readback and therefore no per-step fence)
                              // [1] steps since the last refresh
                              // [2] refresh count, [3] displacement-triggered count
  int dual_on = -1;           // -1 = not yet probed, 0 = off, 1 = on
  int dual_every = 10;
  int dual_stats = 0;
  int dual_check = 0;         // >0: every N steps, count in-cutoff pairs in both
                              // lists and require the counts to be equal
  int dual_nall = 0;
  double dual_skin = 1.8;
  double inner_cutsq = 0.0;
  bigint dual_built_step = -1;

  // FULL and DUAL are template parameters, not members: the default
  // (half list, no dual cutoff) instantiation must generate exactly the code
  // it did before either feature existed.  Carrying them as runtime members
  // cost the default path ~7% on A100 -- three extra views in the functor plus
  // a branch in the inner loop is enough to move the register allocation.
  template<bool FULL, bool DUAL> void flat_compute();
  template<class FunctorT> void fill_common(FunctorT &ff);
  template<int CI> void union_compute();
  template<int CI> void union_build();
  void dual_setup(bool &force_refresh);
  void dual_decide(bool force_refresh);
  void dual_verify();
};

}

#endif
#endif

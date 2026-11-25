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
#include "kokkos_base.h"

#include <vector>

namespace LAMMPS_NS {

// packed j-side gather record; see the block comment in the .cpp
struct alignas(16) LJCL2XQ { KK_FLOAT x, y, z, q; };

// cluster-pair list entry (8 i x 4 j tile), see the cluster block comment in the .cpp
//   cj    j-cluster index (4 consecutive cluster slots)
//   mask  bit kj*8+ki set iff slot pair (ki,kj) is in the list
//   sidx  index into the special-code array, or -1 if the tile has no special pair
//   shift periodic image of the j-cluster, (sx+1) + 3*(sy+1) + 9*(sz+1)
struct alignas(16) LJCL2CTile { int cj; unsigned int mask; int sidx; int shift; };
// axis-aligned bounding box of a cluster (at list-build positions)
struct LJCL2BB { float lo[3], hi[3]; };

template<class DeviceType>
class PairLJCutCoulLong2Kokkos : public PairLJCutCoulLongKokkos<DeviceType>, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairLJCutCoulLong2Kokkos(class LAMMPS *);
  ~PairLJCutCoulLong2Kokkos() override;

  void compute(int, int) override;
  void init_style() override;
  void finish_async_compute() override;
  bool async_compute_pending() const override { return cl_pending && cl_kfirst; }

  // cluster-pair path; public only because nvcc forbids extended device
  // lambdas inside private/protected member functions
  void cl_build();
  void cl_step(int eflag, int vflag);
  template<bool RAT, bool EV> void cl_force();
  template<bool RAT, bool EV, bool UNICUT> void cl_force_launch();

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

  // ---- warp-aggregated inner-list append, LMP_LJCL2_BALLOT=1.
  // Measured (kokkos_neigh.md 10.7): a refresh step costs 5.072 ms against
  // 3.000 for a plain flat step.  The master walk accounts for the 3.000, so
  // the remaining 2.07 ms is the append itself -- ~100M contended
  // atomic_fetch_add on d_innum(i), where every vector lane of an atom hits the
  // SAME counter and the hardware serializes them.
  //
  // d_innum(i) is only ever touched by the one group of vector lanes that owns
  // atom i, so that group can agree among itself first and touch memory once:
  // __match_any_sync supplies the group mask, __ballot_sync the keeper mask,
  // then a single lane reserves the group's whole block with one atomic and
  // broadcasts the base.  Atomics per row fall from one per keeper (~360 at
  // inner skin 2.5) to one per iteration that keeps anything (~95), and no two
  // lanes of a group contend for the same address any more.
  //
  // The base deliberately comes from an atomic rather than from a running count
  // carried in registers: __activemask() carries no convergence guarantee, and a
  // register-carried count would silently drop pairs if a group were ever split
  // across two masks (see the long comment on the append itself).
  //
  // This is NOT the per-lane segmented append that lost 10% in section 7: the
  // row stays compact (so the tight walk is unchanged and still unrollable), the
  // loop is not split, and the body is not inlined twice, so the register
  // allocation should not move.
  int ballot_on = -1;         // -1 = not yet probed, 0 = off, 1 = on

  // ---- packed j-side gather (experimental, LMP_LJCL2_PACK=1) ----
  // Per neighbor entry the kernel gathers x(j,0..2) (three separate 32-bit
  // loads: a LayoutRight float*[3] has a 12-byte stride, which nvcc cannot
  // vectorize), type(j) and, inside the cutoff, q(j).  A per-step packed copy
  // of (x,y,z,q) as one 16-byte-aligned record plus a byte-wide type array
  // turns that into a single LDG.128 and one byte load.  The pack pass is a
  // fully coalesced ~16 MB/step streaming kernel (~15 us).  See kokkos_neigh.md
  // section 18.2 for the request-rate diagnosis this targets.
  Kokkos::View<LJCL2XQ*, DeviceType> d_xq;
  Kokkos::View<unsigned char*, DeviceType> d_type8;
  int pack_on = -1;           // -1 = not yet probed, 0 = off, 1 = on

  // ---- kernel variant (LMP_LJCL2_VARIANT), see the enum in the .cpp.  0 is the
  // real kernel, 1-3 delete one component each to measure its cost directly
  // (their forces are wrong), 4 is the intrinsic-math candidate.
  int variant = -1;           // -1 = not yet probed

  // ---- warp-compacted kernel (experimental, LMP_LJCL2_COMPACT=1, CUDA only).
  // Test phase and physics phase are split by a per-warp shared-memory queue
  // so the physics runs on fully converged warps; see the block comment in
  // the .cpp.  LMP_LJCL2_SMEM=1 moves the type tables to shared memory (A/B;
  // slower on A100), LMP_LJCL2_APW sets the atoms walked per warp.  Combines
  // with LMP_LJCL2_PACK and LMP_LJCL2_DUAL.
  int compact_on = -1;        // -1 = not yet probed
  int compact_smem = 0;
  int compact_apw = 8;
  Kokkos::View<KK_FLOAT*, DeviceType> d_tab;   // packed per-type-pair table
  std::vector<KK_FLOAT> h_tab;                 // last uploaded copy
  void tab_refresh();
  template<bool FULL, bool RAT, bool SMEM, bool PACK = false, bool DUAL = false>
  void compact_compute();

  // ---- shared-memory tiled kernel (experimental, LMP_LJCL2_LOCAL=1|2, CUDA +
  // single precision, half list).  Groups of local_ni consecutive i atoms get
  // a j-set and 16-bit local row entries at each reneighbor; see the block
  // comment in the .cpp.  LMP_LJCL2_LNI sets the group size, pair/team/size
  // the block size (default 256 threads), LMP_LOCAL_STATS prints the tiling.
  int local_on = -1;          // -1 = not yet probed, 0 = off, 1 = float, 2 = fixed point
  int local_ni = 64;
  int local_vec = 1;          // row entries per lane per test pass (LMP_LJCL2_LVEC)
  int local_jcap = 0;
  int local_maxm = 0;
  int local_ngroup = 0;
  int local_stats = 0;
  bigint local_built_step = -1;
  Kokkos::View<unsigned short**, Kokkos::LayoutRight, DeviceType> d_lnbr;
  Kokkos::View<int*, DeviceType> d_ljset;
  Kokkos::View<int*, DeviceType> d_lm;
  Kokkos::View<int*, DeviceType> d_lstat;
  Kokkos::View<unsigned short**, Kokkos::LayoutRight, DeviceType> d_linr;   // LOCAL+DUAL
  Kokkos::View<int*, DeviceType> d_linnum;
  void local_build();
  void local_compute();

  // FULL, DUAL and PACK are template parameters, not members: the default
  // (half list, no dual cutoff) instantiation must generate exactly the code
  // it did before any of them existed.  Carrying them as runtime members
  // cost the default path ~7% on A100 -- three extra views in the functor plus
  // a branch in the inner loop is enough to move the register allocation.
  // ---- cluster-pair list path (experimental, LMP_LJCL2_CLUSTER=1, CUDA, one MPI
  // rank, orthogonal periodic box).  GROMACS-style: atoms are binned into x-y
  // columns, sorted by z and cut into 8-atom i-clusters / 4-atom j-clusters; the
  // tile list is searched directly with bounding boxes (no flat list is built:
  // the pair's neighbor request is made occasional), periodic images use shift
  // vectors instead of ghost atoms, and an inner list pruned to rc + delta is
  // refreshed on the device whenever an atom moved more than delta/2.
  int cl_on = -1;                     // -1 = not probed
  bigint cl_built_step = -1;
  int cl_ncx = 0, cl_ncy = 0, cl_ncol = 0, cl_nslot = 0, cl_nicl = 0;
  int cl_maxout = 0, cl_maxsc = 0;
  double cl_colw[2] = {0.0, 0.0};
  double cl_delta = 1.0;              // inner-list buffer (LMP_CL_DELTA)
  double cl_colscale = 1.0;           // column width scale (LMP_CL_COLSCALE)
  int cl_stats = 0;                   // LMP_CL_STATS
  int cl_nprune = 0, cl_nsteps = 0;
  Kokkos::View<unsigned int*, DeviceType> d_cl_key;
  Kokkos::View<int*, DeviceType> d_cl_perm, d_cl_colcnt, d_cl_colfirst, d_cl_colslot,
    d_cl_atom, d_cl_icol, d_cl_nout, d_cl_nin, d_cl_ctl;
  Kokkos::View<LJCL2XQ*, DeviceType> d_cl_xq;
  Kokkos::View<float*[4], Kokkos::LayoutRight, DeviceType> d_cl_xp;   // x at last prune
  Kokkos::View<unsigned char*, DeviceType> d_cl_t8;
  Kokkos::View<KK_ACC_FLOAT*[3], Kokkos::LayoutRight, DeviceType> d_cl_f;
  Kokkos::View<LJCL2BB*, DeviceType> d_cl_jbb;
  Kokkos::View<LJCL2CTile**, Kokkos::LayoutRight, DeviceType> d_cl_out, d_cl_in;
  Kokkos::View<unsigned int*[2], Kokkos::LayoutRight, DeviceType> d_cl_sc;
  Kokkos::View<double*, DeviceType> d_cl_ev;

  // LMP_CL_ASYNC: on force-only steps the cluster kernels run on a second
  // stream, overlapping bonded + kspace; finish_async_compute() joins them
  int cl_async = 0;
  int cl_kfirst = 1;          // LMP_CL_KFIRST: ask VerletKokkos to run kspace before bonded
  int cl_persist = 0;         // LMP_CL_PERSIST: force-kernel blocks per SM (0 = one warp per i-cluster)
  int cl_nsm = 0;
  bool cl_pending = false;
  typename DeviceType::execution_space cl_space;      // the cl_stream instance
  typename DeviceType::execution_space cl_launch_space; // instance for cl_force_launch
#if defined(KOKKOS_ENABLE_CUDA)
  cudaStream_t cl_stream = nullptr;
  cudaEvent_t cl_ev_in = nullptr, cl_ev_out = nullptr;
#endif
  void cl_scatter();

  int launch_vector_length() const;
  void pack_refresh();
  template<bool FULL, bool DUAL, bool PACK, int VAR = 0, bool BALLOT = false>
  void flat_compute();
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

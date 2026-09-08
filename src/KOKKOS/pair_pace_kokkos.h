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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pace/kk,PairPACEKokkos<LMPDeviceType>);
PairStyle(pace/kk/device,PairPACEKokkos<LMPDeviceType>);
#ifdef LMP_KOKKOS_GPU
PairStyle(pace/kk/host,PairPACEKokkos<LMPHostType>);
#else
PairStyle(pace/kk/host,PairPACEKokkos<LMPDeviceType>);
#endif
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_PACE_KOKKOS_H
#define LMP_PAIR_PACE_KOKKOS_H

#include "pair_pace.h"
#include "kokkos_type.h"
#include "pair_kokkos.h"

class SplineInterpolator;

namespace LAMMPS_NS {

template<class DeviceType>
class PairPACEKokkos : public PairPACE {
 public:
  struct TagPairPACEComputeNeigh{};
  struct TagPairPACEComputeAi{};
  struct TagPairPACEComputeRho{};
  struct TagPairPACEComputeFS{};
  struct TagPairPACEComputeWeights{};
  struct TagPairPACEComputeDerivative{};

  // CPU backend variants: one atom per thread, neighbours and basis functions
  // looped over serially inside the kernel.  That removes every atomic and
  // makes the per-atom accumulators the innermost working set.
  struct TagPairPACEComputeAiCPU{};
  struct TagPairPACEConjugateAi{};   // host only: builds the full A from A_sph
  struct TagPairPACEComputeRadialCPU{};
  struct TagPairPACEComputeRhoCPU{};   // fused rho -> F(rho) -> weights
  struct TagPairPACEComputeRhoBatchCPU{};   // same, 8 atoms per work item
  struct TagPairPACEComputeDerivativeCPU{};

  template<int NEIGHFLAG, int EVFLAG>
  struct TagPairPACEComputeForce{};

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  using complex = SNAComplex<KK_FLOAT>;

  PairPACEKokkos(class LAMMPS *);
  ~PairPACEKokkos() override;

  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeNeigh>::member_type& team) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeAi,const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeAi>::member_type& team) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeRho,const int& iter) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeFS,const int& ii) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeWeights,const int& iter) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeAiCPU,const int& ii) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEConjugateAi,const int& ii) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeRadialCPU,const int& ii) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeRhoCPU,const int& ii) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeRhoBatchCPU,const int& ib) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeDerivativeCPU,const int& ii) const;


// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeDerivative,const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeDerivative>::member_type& team) const;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeForce<NEIGHFLAG,EVFLAG>,const int& ii) const;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeForce<NEIGHFLAG,EVFLAG>,const int& ii, EV_FLOAT&) const;

 protected:
  int inum, maxneigh, chunk_size, chunk_offset, idx_ms_combs_max, idx_sph_max;

  // compile-time host/device switch, so the branches below are "if constexpr"
  // and the unused kernel set is never instantiated (cf. pair_snap_kokkos.h)
  static constexpr int host_flag =
      (ExecutionSpaceFromDevice<DeviceType>::space == LAMMPS_NS::HostKK);

  // set by init_style() when a CPU backend has to defer to the non-accelerated
  // base class because there are no KOKKOS kernels for the requested case
  int host_fallback;

  // team scratch memory level used by the ComputeNeigh short neighbor list
  // build, chosen from the base class neigh_scratch_request (the "neigh"
  // pair_style keyword) and what the device can actually provide
  int neigh_scratch_level;    // level actually used by ComputeNeigh (0 or 1)
  int neigh_scratch_warned;   // whether the auto-fallback warning was printed

  int neigh_scratch_level_select(int scratch_size, int max_level0);

  int eflag, vflag;

  int neighflag, max_ndensity;
  int nelements, lmax, nradmax, nradbase;

  // ------------------------------------------------------------------
  // GPU performance tuning constants (per-architecture).
  //
  // These set the team size used to launch each TeamPolicy kernel. They
  // replace a single hard-coded team size (formerly 32 for every GPU
  // kernel) so that occupancy can be tuned per kernel and per backend.
  //
  // They are selected by host_flag, not by the build configuration: in a GPU
  // build the host instantiation (pace/kk/host) is a CPU backend and wants
  // team size 1. Keying them off KOKKOS_ENABLE_CUDA and friends would give it
  // team size 32, which for ComputeNeigh also means requesting 32x the team
  // scratch it actually needs.
  //
  // NOTE: the GPU values intentionally reproduce the previous behaviour
  // (team size 32 on all GPU kernels). They are the hook for the empirical
  // per-architecture tuning sweep; do not change them without benchmarking.
  // ------------------------------------------------------------------
  static constexpr int team_size_compute_neigh = host_flag ? 1 : 32;
  static constexpr int team_size_compute_ai = host_flag ? 1 : 32;
  static constexpr int team_size_compute_derivative = host_flag ? 1 : 32;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;

  typedef Kokkos::DualView<KK_FLOAT**, DeviceType> tdual_fparams;
  tdual_fparams k_cutsq, k_scale;
  typedef Kokkos::View<KK_FLOAT**, DeviceType> t_fparams;
  t_fparams d_cutsq, d_scale;
  t_fparams d_cut_in, d_dcut_in; // inner cutoff

  typename AT::t_int_1d d_map;

  int need_dup;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;

  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  friend void pair_virial_fdotr_compute<PairPACEKokkos>(PairPACEKokkos*);

  void grow(int, int);
  void copy_pertype();
  void copy_splines();
  void copy_tilde();
  void allocate() override;
  void precompute_harmonics();
  double memory_usage() override;

  template<int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &fx, const KK_FLOAT &fy, const KK_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void cutoff_func_poly(const KK_FLOAT, const KK_FLOAT, const KK_FLOAT, KK_FLOAT &, KK_FLOAT &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void Fexp(const KK_FLOAT, const KK_FLOAT, KK_FLOAT &, KK_FLOAT &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void FexpShiftedScaled(const KK_FLOAT, const KK_FLOAT, KK_FLOAT &, KK_FLOAT &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void inner_cutoff(const KK_FLOAT, const KK_FLOAT, const KK_FLOAT, KK_FLOAT &, KK_FLOAT &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void FS_values_and_derivatives(const int, KK_FLOAT&, const int) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void evaluate_splines(const int, const int, KK_FLOAT, int, int, int, int) const;

  // Shared inner radial loop for ComputeDerivative. Accumulates the gradient
  // contribution of a single (l, m) spherical-harmonic channel into f_ji for
  // all radial functions n. wscale folds in the factor-of-2 used for m > 0.
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void compute_derivative_radial(const int ii, const int jj, const int mu_j,
      const int idx_sph, const int l, const complex &ylm, const complex (&dylm)[3],
      const KK_FLOAT rinv, const KK_FLOAT (&r_hat)[3], const KK_FLOAT wscale,
      KK_ACC_FLOAT (&f_ji)[3]) const;

  // Read A(l,m) for atom-slot ii, element mu, radial index n from the half-basis
  // A_sph using conjugate symmetry A(l,-p) = (-1)^p * conj(A(l,p)). Shared by
  // ComputeRho (product) and ComputeWeights (leave-one-out product recompute).
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  complex read_A(const int ii, const int mu, const int l, const int m, const int n) const;
  // shared body of ComputeAi.  NEED_ATOMICS is true when several threads may
  // accumulate into the same atom (device backends), false when one thread
  // owns the atom (CPU backends), following the sna_kokkos pattern.
  template<bool NEED_ATOMICS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void compute_ai_one(const int, const int) const;

  // Base pointers and strides for one atom, hoisted out of the basis-function
  // loops so the multi-dimensional View subscripts are not recomputed for
  // every access.  Host only, where the layout is LayoutRight.
  struct BasisPtrs {
    const complex *A;
    const KK_FLOAT *A_rank1;
    complex *dB;
    complex *w;
    KK_FLOAT *w_rank1;
    KK_FLOAT *rho;
    const KK_FLOAT *dF;
    const int *mus, *ns, *ls, *ms, *idx_funcs, *rank, *idx_sph;
    const int *A_off, *w_off, *wm_off, *r1_off, *dB_off;
    const KK_FLOAT *ctildes;
    int A_l, A_n, w_l, w_n, rankmax, ndensitymax;
  };

  KOKKOS_INLINE_FUNCTION
  void set_basis_ptrs(BasisPtrs &, const int, const int) const;

  // CPU-only bodies for one (atom, ms-combination): the rank-length products
  // live on the stack instead of in chunk-sized global arrays.  NDENSITY
  // fixes the density count at compile time (0 = runtime).
  template<int NDENSITY>
  KOKKOS_INLINE_FUNCTION
  void rho_fs_weights_cpu(const int) const;

  template<int NDENSITY>
  KOKKOS_INLINE_FUNCTION
  void rho_fs_weights_batch_cpu(const int, const int, const int) const;

  template<int NDENSITY>
  KOKKOS_INLINE_FUNCTION
  void compute_rho_one_cpu(const BasisPtrs &, const int, const int) const;

  template<int NDENSITY>
  KOKKOS_INLINE_FUNCTION
  void compute_weights_one_cpu(const BasisPtrs &, const int, const int) const;

  template<int NDENSITY>
  KOKKOS_INLINE_FUNCTION
  void rho_one_rank1_cpu(const BasisPtrs &, const int, const int) const;

  template<int NDENSITY>
  KOKKOS_INLINE_FUNCTION
  void weights_one_rank1_cpu(const BasisPtrs &, const int, const int) const;

  KOKKOS_INLINE_FUNCTION
  void compute_fs_one(const int) const;

  KOKKOS_INLINE_FUNCTION
  void compute_derivative_one(const int, const int) const;

  // upper bound on the ACE correlation order, for the stack temporaries above
  static constexpr int MAX_RANK_CPU = 16;

  // flat one-atom-per-thread policy used by the CPU kernels.  A dynamic
  // schedule matters because the per-atom cost follows the neighbor count,
  // which is ragged (cf. snap_get_policy in pair_snap_kokkos.h).
  template<class TagStyle>
  auto host_atom_policy() const {
    return Kokkos::RangePolicy<DeviceType, Kokkos::Schedule<Kokkos::Dynamic>,
                               TagStyle>(0, chunk_size);
  }

  template<class TagStyle>
  void check_team_size_for(int, int&, int);

  template<class TagStyle>
  void check_team_size_reduce(int, int&, int);

  // Utility routine which wraps computing per-team scratch size requirements for
  // ComputeNeigh, ComputeUi, and ComputeFusedDeidrj
  template <typename scratch_type>
  int scratch_size_helper(int values_per_team);

  typedef Kokkos::View<int*, DeviceType> t_ace_1i;
  typedef Kokkos::View<int**, DeviceType> t_ace_2i;
  typedef Kokkos::View<int**, Kokkos::LayoutRight, DeviceType> t_ace_2i_lr;
  typedef Kokkos::View<int***, DeviceType> t_ace_3i;
  typedef Kokkos::View<int***, Kokkos::LayoutRight, DeviceType> t_ace_3i_lr;
  typedef Kokkos::View<int****, DeviceType> t_ace_4i;
  typedef Kokkos::View<KK_FLOAT*, DeviceType> t_ace_1d;
  typedef Kokkos::View<KK_FLOAT**, DeviceType> t_ace_2d;
  typedef Kokkos::View<KK_FLOAT**, Kokkos::LayoutRight, DeviceType> t_ace_2d_lr;
  typedef Kokkos::View<KK_FLOAT*[3], DeviceType> t_ace_2d3;
  typedef Kokkos::View<KK_FLOAT***, DeviceType> t_ace_3d;
  typedef Kokkos::View<KK_FLOAT**[3], DeviceType> t_ace_3d3;
  typedef Kokkos::View<KK_FLOAT**[4], DeviceType> t_ace_3d4;
  typedef Kokkos::View<KK_FLOAT**[4], Kokkos::LayoutRight, DeviceType> t_ace_3d4_lr;
  typedef Kokkos::View<KK_FLOAT****, DeviceType> t_ace_4d;
  typedef Kokkos::View<complex*, DeviceType> t_ace_1c;
  typedef Kokkos::View<complex**, DeviceType> t_ace_2c;
  typedef Kokkos::View<complex***, DeviceType> t_ace_3c;
  typedef Kokkos::View<complex**[3], DeviceType> t_ace_3c3;
  typedef Kokkos::View<complex****, DeviceType> t_ace_4c;
  typedef Kokkos::View<complex***[3], DeviceType> t_ace_4c3;

  typedef typename Kokkos::View<KK_FLOAT*, DeviceType>::host_mirror_type th_ace_1d;

  t_ace_3d A_rank1;


  t_ace_3d weights_rank1;
  // Spherical-harmonic weights, stored as separate real/imaginary arrays
  // (rather than an interleaved complex array) so that the atomic
  // accumulation in ComputeWeights is coalesced across the (innermost) atom
  // index on GPUs. Mirrors the ulisttot_re/ulisttot_im layout in Kokkos SNAP.
  t_ace_4d weights_re;
  t_ace_4d weights_im;

  t_ace_1d e_atom;
  t_ace_2d rhos;
  t_ace_2d dF_drho;


  // hard-core repulsion
  t_ace_1d rho_core;
  t_ace_2d cr;
  t_ace_2d dcr;
  t_ace_1d dF_drho_core;
  t_ace_1d dF_dfcut;
  t_ace_1d d_corerep;
  th_ace_1d h_corerep;

  // radial functions
  t_ace_3d fr;
  t_ace_3d dfr;
  t_ace_3d gr;
  t_ace_3d dgr;
  t_ace_3d d_values;
  t_ace_3d d_derivatives;

  // Spherical Harmonics

  void pre_compute_harmonics(int);

  // Spherical-harmonic basis A, stored as separate real/imaginary arrays
  // (rather than an interleaved complex array) so that the atomic
  // accumulation over neighbors in ComputeAi is coalesced across the
  // (innermost) atom index on GPUs. Mirrors Kokkos SNAP's ulisttot_re/im.
  // ------------------------------------------------------------------
  // Host-only ("dual layout") arrays.
  //
  // The CPU kernels keep the original interleaved-complex layout, which is
  // what a single thread walking one atom wants: re/im land in the same
  // cache line. The device kernels use the split re/im arrays below, which
  // is what a warp scattering across atoms wants (coalescing). Only one set
  // is ever allocated -- see grow() -- so the duplication is in the
  // declarations, not in memory.
  // ------------------------------------------------------------------
  t_ace_4c A;           // full (l,m) A basis, host only
  t_ace_4c A_sph;       // interleaved half-basis A, host only
  t_ace_4c weights;     // interleaved weights, host only
  t_ace_3c dB_flatten;  // stored leave-one-out products, host only

  t_ace_4d A_sph_re;
  t_ace_4d A_sph_im;
  t_ace_1d d_idx_sph;
  t_ace_1i d_idx_sph_cpu;   // int copy, sentinel -1 remapped to a trash row
  t_ace_1d alm;
  t_ace_1d blm;
  t_ace_1d cl;
  t_ace_1d dl;

  // short neigh list
  t_ace_1i d_ncount;
  t_ace_2d d_mu;
  t_ace_2d d_rnorms;
  t_ace_3d3 d_rhats;
  t_ace_2i d_nearest;

  // for ZBL core-rep implementation
  t_ace_1d  d_d_min; // [i] -> min-d for atom ii, d=d = r - (cut_in(mu_i, mu_j) - dcut_in(mu_i, mu_j))
  t_ace_1i  d_jj_min; // [i] -> jj-index of nearest neigh (by r-(cut_in-dcut_in) criterion)
  bool is_zbl;

  // per-type
  t_ace_1i d_ndensity;
  t_ace_1i d_npoti;
  t_ace_1d d_rho_core_cutoff;
  t_ace_1d d_drho_core_cutoff;
  t_ace_1d d_E0vals;
  t_ace_2d_lr d_wpre;
  t_ace_2d_lr d_mexp;

  // tilde
  t_ace_1i d_idx_ms_combs_count;
  t_ace_1i d_nms_rank1;   // number of rank-1 ms-combinations, which come first
  t_ace_2i_lr d_rank;
  t_ace_2i_lr d_num_ms_combs;
  t_ace_2i_lr d_idx_funcs;
  t_ace_3i_lr d_mus;
  t_ace_3i_lr d_ns;
  t_ace_3i_lr d_ls;
  t_ace_3i_lr d_ms_combs;
  t_ace_3d d_ctildes;

  // CPU only: the flattened offset of each (ms-combination, rank) term into an
  // atom's A and weights blocks, precomputed once so the inner loops do a
  // single load instead of walking mus/ns/ls/ms_combs and rebuilding the
  // subscript.  d_wm_off carries the (l,-m) offset with the sign of the
  // half-basis factor packed into its low bit.
  t_ace_3i_lr d_A_off, d_w_off, d_wm_off;
  t_ace_2i_lr d_r1_off;
  // prefix offsets into the rank-compacted dB store: entry idx begins at
  // d_dB_off(mu,idx).  Compaction halves the dB stream for typical bases,
  // where most ms-combinations have rank ~3 but the padded stride is rankmax.
  t_ace_2i_lr d_dB_off;
  int dB_total_max;
  // whether the basis fits the stack buffers of the batched CPU kernels
  bool use_batched_cpu = false;

  void build_cpu_offset_tables();

  t_ace_3d3 f_ij;

  void deallocate_views_of_views();

 public:
  struct SplineInterpolatorKokkos {
    int ntot, nlut, num_of_functions;
    KK_FLOAT cutoff, deltaSplineBins, invrscalelookup, rscalelookup;

    t_ace_3d4_lr lookupTable;

    void operator=(const SplineInterpolator &spline);

    void deallocate() {
      lookupTable = t_ace_3d4_lr();
    }

    KK_FLOAT memory_usage() {
      return lookupTable.span() * sizeof(typename decltype(lookupTable)::value_type);
    }

// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    void calcSplines(const int ii, const int jj, const KK_FLOAT r, const t_ace_3d &d_values, const t_ace_3d &d_derivatives) const;
  };

  Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType> k_splines_gk;
  Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType> k_splines_rnl;
  Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType> k_splines_hc;

};
}    // namespace LAMMPS_NS

#endif
#endif

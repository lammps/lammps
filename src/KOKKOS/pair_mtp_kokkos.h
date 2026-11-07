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

//
// Contributing author, Richard Meng, Queen's University at Kingston, 21.01.24, contact@richardzjm.com
//

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mtp/kk,PairMTPKokkos<LMPDeviceType>);
PairStyle(mtp/kk/device,PairMTPKokkos<LMPDeviceType>);
PairStyle(mtp/kk/host,PairMTPKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_MTP_KOKKOS_H
#define LMP_PAIR_MTP_KOKKOS_H

#include "pair_mtp.h"

#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {

template <class DeviceType> class PairMTPKokkos : public PairMTP {
 public:
  // Structs for kernels
  struct TagPairMTPComputeAlphaBasic {};
  struct TagPairMTPComputeAlphaTimes {};
  struct TagPairMTPComputeNbhDers {};
  struct TagPairMTPComputeNbhDersLong {};
  template <int NEIGHFLAG, int EVFLAG> struct TagPairMTPComputeForce {};

  static constexpr int ATOM_TILE_SIZE = 32;
  static constexpr int REVERSE_LONG_THRESHOLD = 128;
  static constexpr int MAX_TEAM_SIZE_BASIC = 64;
  static constexpr int MAX_TEAM_SIZE_GRAPH = 4;
  static constexpr int NEIGH_CAPACITY_ALIGN = 32;
  static constexpr int GRAPH_PARTITION_MAX = 16;
  static constexpr int GRAPH_PARTITION_TARGET = 2048;
  // Kokkos caps a team parallel_reduce grid at this many blocks and strides the
  // rest; parallel_for does not.
  static constexpr int FORCE_MAX_BLOCKS = 32768;

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairMTPKokkos(class LAMMPS *);
  ~PairMTPKokkos() override;

  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void prepare_waves();    //Precalculates node waves and rule lists

  // ========== Kokkos kernels ==========

  //Utility routines
  template <typename scratch_type>
  int scratch_size_helper(int values_per_team);    // Helps calcs scratch size for calcalphabasic

  template <int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION void v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
                                          const KK_FLOAT &fx, const KK_FLOAT &fy,
                                          const KK_FLOAT &fz, const KK_FLOAT &delx,
                                          const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  // ---------- MTP routines (in order of execution) ----------
  // Kernels for computation
  KOKKOS_INLINE_FUNCTION
  void
  operator()(TagPairMTPComputeAlphaBasic,
             const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaBasic>::member_type
                 &team) const;

  KOKKOS_INLINE_FUNCTION
  void
  operator()(TagPairMTPComputeAlphaTimes,
             const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaTimes>::member_type
                 &team) const;

  KOKKOS_INLINE_FUNCTION
  void
  operator()(TagPairMTPComputeNbhDers,
             const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeNbhDers>::member_type
                 &team) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(
      TagPairMTPComputeNbhDersLong,
      const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeNbhDersLong>::member_type
          &team) const;

  template <int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION void
  operator()(const TagPairMTPComputeForce<NEIGHFLAG, EVFLAG> &,
             const typename Kokkos::TeamPolicy<
                 DeviceType, TagPairMTPComputeForce<NEIGHFLAG, EVFLAG>>::member_type &team,
             EV_FLOAT &ev) const;

  template <int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION void
  operator()(const TagPairMTPComputeForce<NEIGHFLAG, 0> &,
             const typename Kokkos::TeamPolicy<
                 DeviceType, TagPairMTPComputeForce<NEIGHFLAG, 0>>::member_type &team) const;

 protected:
  template <int NEIGHFLAG, int EVFLAG>
  EV_FLOAT compute_force(const typename DeviceType::execution_space &, int, int, int);

  // Needed to process the computation in batches to avoid running out of VRAM.
  // input_chunk_size (the chunksize keyword) is parse on parent class
  int chunk_size, chunk_offset;
  int inum, max_valid_neighs, num_waves;
  int wave_begin, wave_end, node_partitions;
  int host_flag, neighflag;
  int eflag, vflag;    // Energy and virial flag
  double inv_cutoff_range, cutoff_sum, radial_mult;
  // Occupancy queries are launch-geometry independent, so they are resolved once.
  int ts_basic, ts_times, ts_nbh, ts_nbh_long, ts_force[2][2];

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
  typename AT::t_int_1d d_map;

  // ---------- Device Arrays ----------
  // Alphas indicies
  Kokkos::View<int **, DeviceType> d_alpha_index_basic;    // For constructing the basic alphas.
  Kokkos::View<int **, DeviceType> d_alpha_index_times;    // For combining alphas

  Kokkos::View<int *, Kokkos::HostSpace> h_waves;    // Node wave offsets
  Kokkos::View<int *, DeviceType> d_wave_nodes;
  Kokkos::View<int *, DeviceType> d_forward_offsets, d_forward_rules;
  Kokkos::View<int *, DeviceType> d_reverse_offsets, d_reverse_split;
  Kokkos::View<int *, Kokkos::HostSpace> h_long_waves;
  Kokkos::View<int *, DeviceType> d_long_nodes;
  Kokkos::View<int *[3], Kokkos::LayoutRight, DeviceType> d_reverse_terms;
  Kokkos::View<int *, DeviceType> d_alpha_moment_mapping;    // Maps alphas to the basis functions.

  // The learned coefficients.
  Kokkos::View<KK_FLOAT *, DeviceType> d_radial_basis_coeffs;    // The radial components.
  Kokkos::View<KK_FLOAT *, DeviceType> d_species_coeffs;         // The species-based constants
  Kokkos::View<KK_FLOAT *, DeviceType> d_linear_coeffs;          // Basis coeffs
  Kokkos::View<KK_FLOAT *, DeviceType> d_moment_coeffs;          // Scalar seeds

  // Global working buffers.
  Kokkos::View<int **, DeviceType> d_valid_neighs;
  Kokkos::View<int *, DeviceType> d_num_valid_neighs;

  // Chunk-local radials
  Kokkos::View<KK_FLOAT ***, Kokkos::LayoutLeft, DeviceType> d_radial_vals;
  Kokkos::View<KK_FLOAT ***, Kokkos::LayoutLeft, DeviceType> d_radial_ders;
  Kokkos::View<KK_FLOAT **, Kokkos::LayoutLeft, DeviceType> d_inv_dist;

  // Indexed by atom tile, moment, and atom lane.
  Kokkos::View<KK_FLOAT **[ATOM_TILE_SIZE], Kokkos::LayoutRight, DeviceType> d_moment_tensor_vals;
  Kokkos::View<KK_FLOAT **[ATOM_TILE_SIZE], Kokkos::LayoutRight, DeviceType>
      d_nbh_energy_ders_wrt_moments;

  // Typedefs for shared memory
  typedef Kokkos::View<KK_FLOAT *, typename DeviceType::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      shared_kk_float_1d;    // Used for basic adjoints
  typedef Kokkos::View<KK_FLOAT **[3], typename DeviceType::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      shared_kk_float_3d;    // Used for coord powers
  typedef Kokkos::View<KK_FLOAT **, typename DeviceType::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      shared_kk_float_2d;    // Used for radial basis vals, ders, and dist powers

  int need_dup;

  // ---------- Define the forces, per-atom energy, and virials----------
  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template <typename DataType, typename Layout>
  using DupScatterView =
      KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template <typename DataType, typename Layout>
  using NonDupScatterView =
      KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT *[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT *[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_ACC_FLOAT *[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT *[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  friend void pair_virial_fdotr_compute<PairMTPKokkos>(PairMTPKokkos *);
};

}    // namespace LAMMPS_NS

#endif
#endif

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
// Contributing author, Richard Meng, Queen's University at Kingston, 10.02.25, contact@richardzjm.com
//

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mtp/extrapolation/kk,PairMTPExtrapolationKokkos<LMPDeviceType>);
PairStyle(mtp/extrapolation/kk/device,PairMTPExtrapolationKokkos<LMPDeviceType>);
PairStyle(mtp/extrapolation/kk/host,PairMTPExtrapolationKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_MTP_EXTRAPOLATION_KOKKOS_H
#define LMP_PAIR_MTP_EXTRAPOLATION_KOKKOS_H

#include "pair_mtp_extrapolation.h"

#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {

template <class DeviceType> class PairMTPExtrapolationKokkos : public PairMTPExtrapolation {
 public:
  // Structs for kernels
  struct TagPairMTPComputeAlphaBasic {};
  struct TagPairMTPComputeAlphaTimes {};
  struct TagPairMTPComputeNbhDers {};
  struct TagPairMTPComputeNbhDersLong {};
  struct TagPairMTPReduceCoeffDers {};
  struct TagPairMTPCombineCoeffDers {};
  template <int NEIGHFLAG, int EVFLAG> struct TagPairMTPComputeForce {};

  static constexpr int ATOM_TILE_SIZE = 32;
  static constexpr int REVERSE_LONG_THRESHOLD = 128;

  // Kokkos caps a team parallel_reduce grid at this many blocks and strides the
  // rest; parallel_for does not.
  static constexpr int FORCE_MAX_BLOCKS = 32768;
  static constexpr int COEFF_REDUCE_BLOCK_SIZE = 1024;

  static constexpr int MAX_TEAM_SIZE_BASIC = 64;
  static constexpr int MAX_TEAM_SIZE_GRAPH = 4;

  static constexpr int NEIGH_CAPACITY_ALIGN = 32;

  static constexpr int GRAPH_PARTITION_MAX = 16;
  static constexpr int GRAPH_PARTITION_TARGET = 2048;

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairMTPExtrapolationKokkos(class LAMMPS *);
  ~PairMTPExtrapolationKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void prepare_waves();    //Precalculates node waves and rule lists

  // ========== Kokkos kernels ==========

  //Utility routines
  template <typename scratch_type> int scratch_size_helper(int values_per_team);

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
  KOKKOS_INLINE_FUNCTION
  void
  operator()(TagPairMTPReduceCoeffDers,
             const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPReduceCoeffDers>::member_type
                 &team) const;
  KOKKOS_INLINE_FUNCTION
  void
  operator()(TagPairMTPCombineCoeffDers,
             const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPCombineCoeffDers>::member_type
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
  void evaluate_grades() override;

  template <int NEIGHFLAG, int EVFLAG>
  EV_FLOAT compute_force(const typename DeviceType::execution_space &, int, int, int);

  // Needed to process the computation in batches to avoid running out of VRAM.
  // input_chunk_size (the chunksize keyword) is parse on parent class
  int chunk_size, chunk_offset;
  int inum, max_valid_neighs, num_waves;
  int wave_begin, wave_end, node_partitions;
  int host_flag, neighflag;
  int eflag, vflag;
  bool calculate_grade_this_step;
  double inv_cutoff_range, cutoff_sum, radial_mult;
  int ts_basic, ts_times, ts_nbh, ts_nbh_long, ts_force[2][2];
  int ts_reduce, ts_nbh_grade, ts_cfg_grade;
  int coeff_reduce_blocks;
  int cached_force_team_scratch_size;

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
  Kokkos::View<int **, DeviceType> d_alpha_index_basic;
  Kokkos::View<int **, DeviceType> d_alpha_index_times;
  Kokkos::View<int *, Kokkos::HostSpace> h_waves;
  Kokkos::View<int *, DeviceType> d_wave_nodes;
  Kokkos::View<int *, DeviceType> d_forward_offsets, d_forward_rules;
  Kokkos::View<int *, DeviceType> d_reverse_offsets, d_reverse_split;
  Kokkos::View<int *, Kokkos::HostSpace> h_long_waves;
  Kokkos::View<int *, DeviceType> d_long_nodes;
  Kokkos::View<int *[3], Kokkos::LayoutRight, DeviceType> d_reverse_terms;
  Kokkos::View<int *, DeviceType> d_alpha_moment_mapping;

  // The learned coefficients.
  Kokkos::View<KK_FLOAT *, DeviceType> d_radial_basis_coeffs;
  Kokkos::View<KK_FLOAT *, DeviceType> d_species_coeffs;
  Kokkos::View<KK_FLOAT *, DeviceType> d_linear_coeffs;
  Kokkos::View<KK_FLOAT *, DeviceType> d_moment_coeffs;

  // Global working buffers.
  Kokkos::View<int **, DeviceType> d_valid_neighs;
  Kokkos::View<int *, DeviceType> d_num_valid_neighs;
  Kokkos::View<KK_FLOAT ***, Kokkos::LayoutLeft, DeviceType> d_radial_vals;
  Kokkos::View<KK_FLOAT ***, Kokkos::LayoutLeft, DeviceType> d_radial_ders;
  Kokkos::View<KK_FLOAT **, Kokkos::LayoutLeft, DeviceType> d_inv_dist;
  // Indexed by atom tile, moment, and atom lane.
  Kokkos::View<KK_FLOAT **[ATOM_TILE_SIZE], Kokkos::LayoutRight, DeviceType> d_moment_tensor_vals;
  Kokkos::View<KK_FLOAT **[ATOM_TILE_SIZE], Kokkos::LayoutRight, DeviceType>
      d_nbh_energy_ders_wrt_moments;

  // Extrapolation-grade buffers.
  Kokkos::View<KK_FLOAT ***, Kokkos::LayoutLeft, DeviceType> d_radial_basis_cache;
  Kokkos::View<KK_FLOAT **, Kokkos::LayoutLeft, DeviceType> d_local_coeff_ders;
  Kokkos::View<KK_FLOAT **, DeviceType> d_inverse_active_set;    // Active-set inverse matrix
  Kokkos::View<KK_FLOAT *, DeviceType> d_nbh_extrapolation_grades;    // Per-atom grades
  Kokkos::View<KK_FLOAT *, DeviceType> d_energy_ders_wrt_coeffs;      // Configuration derivative
  Kokkos::View<KK_FLOAT **, Kokkos::LayoutLeft, DeviceType> d_coeff_ders_partials;    // Chunk sums

  // Typedefs for shared memory
  typedef Kokkos::View<KK_FLOAT *, typename DeviceType::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      shared_kk_float_1d;
  typedef Kokkos::View<KK_FLOAT **[3], typename DeviceType::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      shared_kk_float_3d;
  typedef Kokkos::View<KK_FLOAT **, typename DeviceType::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      shared_kk_float_2d;

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

  friend void pair_virial_fdotr_compute<PairMTPExtrapolationKokkos>(PairMTPExtrapolationKokkos *);
};

template <class DeviceType> struct ComputeNbhGrades {
  // Computes grades for a tile of neighbourhoods from local coefficient derivatives.
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef KK_FLOAT value_type;
  static constexpr int ATOM_TILE_SIZE = PairMTPExtrapolationKokkos<DeviceType>::ATOM_TILE_SIZE;
  static constexpr int NBH_TILE_SIZE = 8;
  static constexpr int COEFF_TILE_SIZE = 64;
  typedef Kokkos::View<KK_FLOAT **[ATOM_TILE_SIZE], Kokkos::LayoutRight, DeviceType> MomentView;
  typedef Kokkos::View<KK_FLOAT *[NBH_TILE_SIZE], Kokkos::LayoutRight,
                       typename DeviceType::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      shared_kk_float_2d;

  int chunk_size, chunk_offset;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d d_map;
  int species_count, radial_coeff_count, radial_coeff_count_per_pair, alpha_scalar_count,
      coeff_count;
  Kokkos::View<KK_FLOAT **, Kokkos::LayoutLeft, DeviceType> d_local_coeff_ders;
  MomentView d_moment_tensor_vals;
  Kokkos::View<int *, DeviceType> d_alpha_moment_mapping;
  Kokkos::View<KK_FLOAT **, DeviceType> d_inverse_active_set;
  Kokkos::View<KK_FLOAT *, DeviceType> d_nbh_extrapolation_grades;

  ComputeNbhGrades(int chunk_size_, int chunk_offset_, typename AT::t_int_1d_randomread d_ilist_,
                   typename AT::t_int_1d_randomread type_, typename AT::t_int_1d d_map_,
                   int species_count_, int radial_coeff_count_, int radial_coeff_count_per_pair_,
                   int alpha_scalar_count_, int coeff_count_,
                   Kokkos::View<KK_FLOAT **, Kokkos::LayoutLeft, DeviceType> d_local_coeff_ders_,
                   MomentView d_moment_tensor_vals_,
                   Kokkos::View<int *, DeviceType> d_alpha_moment_mapping_,
                   Kokkos::View<KK_FLOAT **, DeviceType> d_inverse_active_set_,
                   Kokkos::View<KK_FLOAT *, DeviceType> d_nbh_extrapolation_grades_) :
      chunk_size(chunk_size_), chunk_offset(chunk_offset_), d_ilist(d_ilist_), type(type_),
      d_map(d_map_), species_count(species_count_), radial_coeff_count(radial_coeff_count_),
      radial_coeff_count_per_pair(radial_coeff_count_per_pair_),
      alpha_scalar_count(alpha_scalar_count_), coeff_count(coeff_count_),
      d_local_coeff_ders(d_local_coeff_ders_), d_moment_tensor_vals(d_moment_tensor_vals_),
      d_alpha_moment_mapping(d_alpha_moment_mapping_), d_inverse_active_set(d_inverse_active_set_),
      d_nbh_extrapolation_grades(d_nbh_extrapolation_grades_)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<DeviceType>::member_type &team,
                  KK_FLOAT &nbh_max_grade) const;
};

template <class DeviceType> struct ComputeCfgGrade {
  // Computes the maximum configuration grade from the accumulated derivatives.
  typedef DeviceType device_type;
  typedef KK_FLOAT value_type;

  const int coeff_count;
  Kokkos::View<KK_FLOAT *, DeviceType> d_energy_ders_wrt_coeffs;
  Kokkos::View<KK_FLOAT **, DeviceType> d_inverse_active_set;

  ComputeCfgGrade(int coeff_count_, Kokkos::View<KK_FLOAT *, DeviceType> d_energy_ders_wrt_coeffs_,
                  Kokkos::View<KK_FLOAT **, DeviceType> d_inverse_active_set_) :
      coeff_count(coeff_count_), d_energy_ders_wrt_coeffs(d_energy_ders_wrt_coeffs_),
      d_inverse_active_set(d_inverse_active_set_)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<DeviceType>::member_type &team,
                  KK_FLOAT &cfg_max_grade) const;
};

}    // namespace LAMMPS_NS

#endif
#endif

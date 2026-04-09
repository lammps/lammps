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

#ifndef LMP_PAIR_MTP_KOKKOS_EXTRAPOLATION_H
#define LMP_PAIR_MTP_KOKKOS_EXTRAPOLATION_H

#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#include "pair_kokkos.h"
#include "pair_mtp_extrapolation.h"

namespace LAMMPS_NS {

template <class DeviceType> class PairMTPExtrapolationKokkos : public PairMTPExtrapolation {
 public:
  // Structs for kernels
  struct TagPairMTPInitMomentValsDers {};
  struct TagPairMTPInitRadJacobian {};
  struct TagPairMTPInitCoeffDers {};
  struct TagPairMTPComputeAlphaBasic {};
  struct TagPairMTPComputeAlphaBasicRad {};
  struct TagPairMTPComputeAlphaTimes {};
  struct TagPairMTPSetScalarNbhDers {};
  struct TagPairMTPComputeNbhDers {};
  struct TagPairMTPReduceCoeffDers {};
  template <int NEIGHFLAG, int EVFLAG> struct TagPairMTPComputeForce {};

  enum { EnabledNeighFlags = HALF | HALFTHREAD };
  enum { COUL_FLAG = 0 };
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairMTPExtrapolationKokkos(class LAMMPS *);
  ~PairMTPExtrapolationKokkos() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void evaluate_grades();
  void prepare_waves();    //Precalculates the waves of alpha times by dependency

  // ========== Kokkos kernels ==========
  //Utility routines
  template <class TagStyle> void check_team_size_for(int, int &, int);

  template <typename scratch_type>
  int scratch_size_helper(int values_per_team);    // Helps calcs scratch sizes

  template <int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION void v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
                                          const KK_FLOAT &fx, const KK_FLOAT &fy,
                                          const KK_FLOAT &fz, const KK_FLOAT &delx,
                                          const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  // ---------- MTP routines (in order of execution) ----------

  //Kernels for initing working views
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMTPInitMomentValsDers, const int &k, const int &ii) const;

  KOKKOS_INLINE_FUNCTION    // Only runs on steps when we sample extrapolation
      void operator()(TagPairMTPInitRadJacobian, const int &ii, const int &k) const;

  KOKKOS_INLINE_FUNCTION    // Only runs on steps when we sample extrapolation
      void operator()(TagPairMTPInitCoeffDers, const int &kk) const;

  // Kernels for computation
  KOKKOS_INLINE_FUNCTION
  void
  operator()(TagPairMTPComputeAlphaBasic,
             const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaBasic>::member_type
                 &team) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(
      TagPairMTPComputeAlphaBasicRad,
      const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaBasicRad>::member_type
          &team) const;

  KOKKOS_INLINE_FUNCTION
  void
  operator()(TagPairMTPComputeAlphaTimes,
             const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaTimes>::member_type
                 &team) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMTPSetScalarNbhDers, const int &k, const int &ii) const;

  KOKKOS_INLINE_FUNCTION
  void
  operator()(TagPairMTPComputeNbhDers,
             const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeNbhDers>::member_type
                 &team) const;

  template <int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION void
  operator()(const TagPairMTPComputeForce<NEIGHFLAG, EVFLAG> &,
             const typename Kokkos::TeamPolicy<
                 DeviceType, TagPairMTPComputeForce<NEIGHFLAG, EVFLAG>>::member_type &team,
             EV_FLOAT &ev) const;

  // Kernels for extrapolation computation
  KOKKOS_INLINE_FUNCTION
  void
  operator()(TagPairMTPReduceCoeffDers,
             const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPReduceCoeffDers>::member_type
                 &team) const;

 protected:
  int input_chunk_size, chunk_size,
      chunk_offset;    // Needed to process the computation in batches to avoid running out of VRAM.

  int inum, max_neighs, max_valid_neighs, num_waves;
  int host_flag, neighflag;

  int eflag, vflag;    // Energy and virial flag

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

  // ---------- Device Arrays  ----------
  // Alphas indicies
  Kokkos::View<int **, DeviceType> d_alpha_index_basic;      // For constructing the basic alphas.
  Kokkos::View<int **, DeviceType> d_alpha_index_times;      // For combining alphas
  Kokkos::View<int *, DeviceType> d_waves;                   // Dependency waves
  Kokkos::View<int *, DeviceType> d_alpha_moment_mapping;    // Maps alphas to the basis functions.

  // The learned coefficients.
  Kokkos::View<KK_FLOAT *, DeviceType> d_radial_basis_coeffs;    // The radial components.
  Kokkos::View<KK_FLOAT *, DeviceType> d_species_coeffs;         // The species-based constants
  Kokkos::View<KK_FLOAT *, DeviceType> d_linear_coeffs;          // Basis coeffs

  // Inverse active set and grades. Only needed in neigbhourhood mode and single MPI configuration mode.
  Kokkos::View<KK_FLOAT **, DeviceType> d_inverse_active_set;
  Kokkos::View<KK_FLOAT *, DeviceType> d_nbh_extrapolation_grades;

  // Source for candidate vector reduction. Only needed in configuration mode.
  // We need two arrays to account for multiple chunks.
  Kokkos::View<KK_FLOAT *, DeviceType> d_energy_ders_wrt_coeffs;
  Kokkos::View<KK_FLOAT *, DeviceType> d_tmp_energy_ders_wrt_coeffs;

  // Global working buffers.
  Kokkos::View<int **, DeviceType> d_valid_neighs;
  Kokkos::View<int *, DeviceType> d_num_valid_neighs;
  Kokkos::View<KK_FLOAT ****, DeviceType> d_moment_jacobian;
  Kokkos::View<KK_FLOAT ***, DeviceType> d_radial_jacobian;
  Kokkos::View<KK_FLOAT **, Kokkos::LayoutRight, DeviceType> d_moment_tensor_vals;
  Kokkos::View<KK_FLOAT **, Kokkos::LayoutRight, DeviceType> d_nbh_energy_ders_wrt_moments;
  Kokkos::View<bool **, DeviceType> d_within_cutoff;

  // Typedefs for shared memory
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

  friend void pair_virial_fdotr_compute<PairMTPExtrapolationKokkos>(PairMTPExtrapolationKokkos *);
};

// ========== Additional functors for grade reductions ==========
template <class DeviceType> struct sComputeNbhGrades {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef KK_FLOAT value_type;

  typedef Kokkos::View<KK_FLOAT *, typename DeviceType::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      shared_kk_float_1d;    // Used for storing coeff derivatives

  const int chunk_size;
  const int chunk_offset;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread type;
  const int species_count;
  const int radial_coeff_count;
  const int alpha_index_basic_count;
  const int radial_coeff_count_per_pair;
  const int alpha_scalar_count;
  const int coeff_count;

  Kokkos::View<KK_FLOAT **, Kokkos::LayoutRight, DeviceType> d_nbh_energy_ders_wrt_moments;
  Kokkos::View<KK_FLOAT ***, DeviceType> d_radial_jacobian;
  Kokkos::View<KK_FLOAT **, Kokkos::LayoutRight, DeviceType> d_moment_tensor_vals;
  Kokkos::View<int *, DeviceType> d_alpha_moment_mapping;
  Kokkos::View<KK_FLOAT **, DeviceType> d_inverse_active_set;
  Kokkos::View<KK_FLOAT *, DeviceType> d_nbh_extrapolation_grades;

  sComputeNbhGrades(
      int chunk_size_, int chunk_offset_, typename AT::t_int_1d_randomread d_ilist_,
      typename AT::t_int_1d_randomread type_, int species_count_, int radial_coeff_count_,
      int alpha_index_basic_count_, int radial_coeff_count_per_pair_, int alpha_scalar_count_,
      int coeff_count_,
      Kokkos::View<KK_FLOAT **, Kokkos::LayoutRight, DeviceType> d_nbh_energy_ders_wrt_moments_,
      Kokkos::View<KK_FLOAT ***, DeviceType> d_radial_jacobian_,
      Kokkos::View<KK_FLOAT **, Kokkos::LayoutRight, DeviceType> d_moment_tensor_vals_,
      Kokkos::View<int *, DeviceType> d_alpha_moment_mapping_,
      Kokkos::View<KK_FLOAT **, DeviceType> d_inverse_active_set_,
      Kokkos::View<KK_FLOAT *, DeviceType> d_nbh_extrapolation_grades_) :
      chunk_size(chunk_size_), chunk_offset(chunk_offset_), d_ilist(d_ilist_), type(type_),
      species_count(species_count_), radial_coeff_count(radial_coeff_count_),
      alpha_index_basic_count(alpha_index_basic_count_),
      radial_coeff_count_per_pair(radial_coeff_count_per_pair_),
      alpha_scalar_count(alpha_scalar_count_), coeff_count(coeff_count_),
      d_nbh_energy_ders_wrt_moments(d_nbh_energy_ders_wrt_moments_),
      d_radial_jacobian(d_radial_jacobian_), d_moment_tensor_vals(d_moment_tensor_vals_),
      d_alpha_moment_mapping(d_alpha_moment_mapping_), d_inverse_active_set(d_inverse_active_set_),
      d_nbh_extrapolation_grades(d_nbh_extrapolation_grades_)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<DeviceType>::member_type &team,
                  KK_FLOAT &nbh_max_grade) const;
};

template <class DeviceType> struct sComputeCfgGrade {
  typedef DeviceType device_type;
  typedef KK_FLOAT value_type;

  const int coeff_count;
  Kokkos::View<KK_FLOAT *, DeviceType> d_energy_ders_wrt_coeffs;
  Kokkos::View<KK_FLOAT **, DeviceType> d_inverse_active_set;

  sComputeCfgGrade(int coeff_count_, Kokkos::View<KK_FLOAT *, DeviceType> d_energy_ders_wrt_coeffs_,
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
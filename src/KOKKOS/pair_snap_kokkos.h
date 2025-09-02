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
PairStyle(snap/kk,PairSNAPKokkosDevice<LMPDeviceType>);
PairStyle(snap/kk/device,PairSNAPKokkosDevice<LMPDeviceType>);
#ifdef LMP_KOKKOS_GPU
PairStyle(snap/kk/host,PairSNAPKokkosHost<LMPHostType>);
#else
PairStyle(snap/kk/host,PairSNAPKokkosDevice<LMPHostType>);
#endif
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_SNAP_KOKKOS_H
#define LMP_PAIR_SNAP_KOKKOS_H

#include "pair_snap.h"
#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {
// pre-declare so sna_kokkos.h can refer to it
template<class DeviceType, typename real_type_, int vector_length_> class PairSNAPKokkos;
};

#include "sna_kokkos.h"

namespace LAMMPS_NS {

// Routines for both the CPU and GPU backend
struct TagPairSNAPPreUi{};
struct TagPairSNAPTransformUi{}; // re-order ulisttot from SoA to AoSoA, zero ylist
template <bool chemsnap> struct TagPairSNAPComputeZi{};
template <bool chemsnap> struct TagPairSNAPComputeBi{};
struct TagPairSNAPComputeBetaLinear{};
struct TagPairSNAPComputeBetaQuadratic{};
template <bool chemsnap> struct TagPairSNAPComputeYi{};
template <bool chemsnap> struct TagPairSNAPComputeYiWithZlist{};
template<int NEIGHFLAG, int EVFLAG>
struct TagPairSNAPComputeForce{};

// GPU backend only
struct TagPairSNAPComputeNeigh{};
struct TagPairSNAPComputeCayleyKlein{};
template <bool chemsnap> struct TagPairSNAPComputeUiSmall{}; // more parallelism, more divergence
template <bool chemsnap> struct TagPairSNAPComputeUiLarge{}; // less parallelism, no divergence
template<int dir> struct TagPairSNAPComputeFusedDeidrjSmall{}; // more parallelism, more divergence
template<int dir> struct TagPairSNAPComputeFusedDeidrjLarge{}; // less parallelism, no divergence
struct TagPairSNAPComputeFusedDeidrjAllSmall{};                // more parallelism, more divergence
struct TagPairSNAPComputeFusedDeidrjAllLarge{};                // less parallelism, no divergence

// CPU backend only
struct TagPairSNAPComputeNeighCPU{};
struct TagPairSNAPComputeUiCPU{};
struct TagPairSNAPComputeDuidrjCPU{};
struct TagPairSNAPComputeDeidrjCPU{};

template<class DeviceType, typename real_type_, int vector_length_>
class PairSNAPKokkos : public PairSNAP {
 public:
  enum {EnabledNeighFlags=FULL|HALF|HALFTHREAD};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  static constexpr int host_flag = (ExecutionSpaceFromDevice<DeviceType>::space == LAMMPS_NS::HostKK);
  static constexpr bool legacy_on_gpu = false; // run the CPU path on the GPU
  static constexpr int vector_length = vector_length_;
  using real_type = real_type_;
  using complex = SNAComplex<real_type>;

  // Static team/tile sizes for device offload

#ifdef KOKKOS_ENABLE_HIP
  static constexpr int team_size_compute_neigh = 2;
  static constexpr int tile_size_compute_ck = 2;
  static constexpr int tile_size_pre_ui = 2;
  static constexpr int base_team_size_compute_ui = 2;
  static constexpr int tile_size_transform_ui = 2;
  static constexpr int tile_size_compute_zi = 2;
  static constexpr int min_blocks_compute_zi = 0; // no minimum bound
  static constexpr int tile_size_compute_bi = 2;
  static constexpr int tile_size_compute_beta = 2;
  static constexpr int tile_size_compute_yi = 2;
  static constexpr int min_blocks_compute_yi = 0; // no minimum bound
  static constexpr int team_size_compute_fused_deidrj = 2;
  static constexpr int team_size_compute_fused_deidrj_all = 1;

  static constexpr int padding_factor = host_flag ? 1 : 4; // extra padding factor
  static constexpr int ui_batch = host_flag ? 1 : 2;
  static constexpr int yi_batch = host_flag ? 1 : 4;
  static constexpr bool use_deidrj_all = true; // whether or not to use the directionally fused deidrj
#elif defined(KOKKOS_ENABLE_SYCL)
  static constexpr int team_size_compute_neigh = 4;
  static constexpr int tile_size_compute_ck = 4;
  static constexpr int tile_size_pre_ui = 8;
  static constexpr int base_team_size_compute_ui = 8;
  static constexpr int tile_size_transform_ui = 8;
  static constexpr int tile_size_compute_zi = 4;
  static constexpr int min_blocks_compute_zi = 0; // no minimum bound
  static constexpr int tile_size_compute_bi = 4;
  static constexpr int tile_size_compute_beta = 8;
  static constexpr int tile_size_compute_yi = 8;
  static constexpr int min_blocks_compute_yi = 0; // no minimum bound
  static constexpr int team_size_compute_fused_deidrj = 4;
  static constexpr int team_size_compute_fused_deidrj_all = 1;

  static constexpr int padding_factor = host_flag ? 1 : 2; // extra padding factor
  static constexpr int ui_batch = 1;
  static constexpr int yi_batch = host_flag ? 1 : 2;
  static constexpr bool use_deidrj_all = false; // whether or not to use the directionally fused deidrj
#else
  static constexpr int team_size_compute_neigh = 4;
  static constexpr int tile_size_compute_ck = 4;
  static constexpr int tile_size_pre_ui = 4;
  static constexpr int base_team_size_compute_ui = sizeof(real_type) == 4 ? 8 : 4;
  static constexpr int tile_size_transform_ui = 4;
  static constexpr int tile_size_compute_zi = 8;
  static constexpr int tile_size_compute_bi = 4;
  static constexpr int tile_size_compute_beta = 4;
  static constexpr int tile_size_compute_yi = 8;
  static constexpr int team_size_compute_fused_deidrj = sizeof(real_type) == 4 ? 4 : 2;
  static constexpr int team_size_compute_fused_deidrj_all = sizeof(real_type) == 4 ? 2 : 1;

  // this empirically reduces perf fluctuations from compiler version to compiler version
  static constexpr int min_blocks_compute_zi = 2;
  static constexpr int min_blocks_compute_yi = 2;

  static constexpr int padding_factor = host_flag ? 1 : 4; // extra padding factor
  static constexpr int ui_batch = host_flag ? 1 : 4;
  static constexpr int yi_batch = host_flag ? 1 : 4;
  static constexpr bool use_deidrj_all = true; // whether or not to use the directionally fused deidrj
#endif

  // Total number of dimensions in the fully fused `ComputeFusedDeidrj`.
  static constexpr int dims = 3;

  // Determine the final batch size for ComputeUi. This convention guarantees that each "team"
  // launched by ComputeUi requests the same amount of scratchpad memory independent of the
  // value of `ui_batch`.
  static constexpr int team_size_compute_ui = base_team_size_compute_ui / ui_batch;
  static_assert(team_size_compute_ui > 0, "ui_batch is too large for team_size_compute_ui");
  static_assert(base_team_size_compute_ui % ui_batch == 0, "ComputeUi batch size must divide into the team size");

  // Check that `yi_batch` evenly divides into `padding_factor`. This guarantees that all data structures are appropriately
  // padded for routines that process `yi_batch` atoms per thread. For the time being, `yi_batch` is re-used across
  // `ComputeYi`, `ComputeZi`, `ComputeYiWithZlist`, and `ComputeBi`. In theory different values could be used across
  // each of these routines, and `padding_factor` would have to be the least common factor of all of these values (or a
  // multiple thereof).
  static_assert((padding_factor / yi_batch) * yi_batch == padding_factor, "yi_batch must divide into padding_factor");

  // Custom MDRangePolicy, Rank3, to reduce verbosity of kernel launches
  // This hides the Kokkos::IndexType<int> and Kokkos::Rank<3...>
  // and reduces the verbosity of the LaunchBound by hiding the explicit
  // multiplication by vector_length
  template <class Device, int num_tiles, class TagPairSNAP, int min_blocks = 0>
  using Snap3DRangePolicy = typename Kokkos::MDRangePolicy<Device, Kokkos::IndexType<int>, Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>, Kokkos::LaunchBounds<vector_length * num_tiles, min_blocks>, TagPairSNAP>;

  // Custom SnapAoSoATeamPolicy to reduce the verbosity of kernel launches
  // This hides the LaunchBounds abstraction by hiding the explicit
  // multiplication by vector length
  template <class Device, int num_teams, class TagPairSNAP>
  using SnapAoSoATeamPolicy = typename Kokkos::TeamPolicy<Device, Kokkos::LaunchBounds<vector_length * num_teams>, TagPairSNAP>;

  // Custom MDRangePolicy, Rank2, on the host, to reduce verbosity of kernel launches. The striding of this launch is intentionally
  // different from the tiled 3D range policy on the device.
  template <class Device, class TagPairSNAP>
  using Snap2DHostRangePolicy = typename Kokkos::MDRangePolicy<Device, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<int>, Kokkos::Rank<2, Kokkos::Iterate::Right, Kokkos::Iterate::Right>, TagPairSNAP>;

  // Custom RangePolicy, Rank2, on the host, to reduce verbosity of kernel launches
  template <class Device, class TagPairSNAP>
  using Snap1DHostRangePolicy = typename Kokkos::RangePolicy<Device, Kokkos::Schedule<Kokkos::Dynamic>, TagPairSNAP>;

  // Helper routine that returns a CPU or a GPU policy as appropriate
  template <class Device, int num_tiles, class TagPairSNAP, int min_blocks = 0>
  auto snap_get_policy(const int& chunk_size_div, const int& second_loop) {
    if constexpr (host_flag) {
      return Snap1DHostRangePolicy<Device, TagPairSNAP>(0, chunk_size_div * vector_length);

      // the 2-d policy is still correct but it has atomics so it's slower on the CPU
      //return Snap2DHostRangePolicy<Device, TagPairSNAP>({0, 0}, {chunk_size_div * vector_length, second_loop});
    } else
      return Snap3DRangePolicy<Device, num_tiles, TagPairSNAP, min_blocks>({0, 0, 0},
                                                                   {vector_length, second_loop, chunk_size_div},
                                                                   {vector_length, num_tiles, 1});
  }

  // Helper routine that dispatches Ui with and without chemsnap
  template <template <bool> class TagPairSNAPComputeUi>
  void snap_dispatch_ui(int n_teams_div, int scratch_size) {
    // make sure we're only passing in types we expect
    static_assert(std::is_same_v<TagPairSNAPComputeUiLarge<false>, TagPairSNAPComputeUi<false>> ||
      std::is_same_v<TagPairSNAPComputeUiSmall<false>, TagPairSNAPComputeUi<false>>);

    std::string name = [&] () -> std::string {
      if constexpr (std::is_same_v<TagPairSNAPComputeUiLarge<false>, TagPairSNAPComputeUi<false>>)
        return (nelements > 1) ? "ComputeUiLargeChemsnap" : "ComputeUiLarge";
      else
        return (nelements > 1) ? "ComputeUiSmallChemsnap" : "ComputeUiSmall";
    }();

    if (nelements > 1) {
      SnapAoSoATeamPolicy<DeviceType, team_size_compute_ui, TagPairSNAPComputeUi<true>>
            policy_ui(n_teams_div, team_size_compute_ui, vector_length);
      policy_ui = policy_ui.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
      Kokkos::parallel_for(name, policy_ui, *this);
    } else {
      SnapAoSoATeamPolicy<DeviceType, team_size_compute_ui, TagPairSNAPComputeUi<false>>
            policy_ui(n_teams_div, team_size_compute_ui, vector_length);
      policy_ui = policy_ui.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
      Kokkos::parallel_for(name, policy_ui, *this);
    }
  }

  // Helper routine that dispatches directional ComputeFusedDeidrj
  template <template <int> class TagPairSNAPComputeFusedDeidrj>
  void snap_dispatch_fused_deidrj(int n_teams_div, int scratch_size) {
    // make sure we're only passing in types we expect
    static_assert(std::is_same_v<TagPairSNAPComputeFusedDeidrjLarge<0>, TagPairSNAPComputeFusedDeidrj<0>> ||
      std::is_same_v<TagPairSNAPComputeFusedDeidrjSmall<0>, TagPairSNAPComputeFusedDeidrj<0>>);

    std::string name = { std::is_same_v<TagPairSNAPComputeFusedDeidrjLarge<0>, TagPairSNAPComputeFusedDeidrj<0>> ? "ComputeFusedDeidrjLarge<0>" : "ComputeFusedDeidrjSmall<0>" };

    // x direction
    SnapAoSoATeamPolicy<DeviceType, team_size_compute_fused_deidrj, TagPairSNAPComputeFusedDeidrj<0> > policy_fused_deidrj_x(n_teams_div,team_size_compute_fused_deidrj,vector_length);
    policy_fused_deidrj_x = policy_fused_deidrj_x.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
    Kokkos::parallel_for(name, policy_fused_deidrj_x, *this);

    // y direction
    name[24] = '1';
    SnapAoSoATeamPolicy<DeviceType, team_size_compute_fused_deidrj, TagPairSNAPComputeFusedDeidrj<1> > policy_fused_deidrj_y(n_teams_div,team_size_compute_fused_deidrj,vector_length);
    policy_fused_deidrj_y = policy_fused_deidrj_y.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
    Kokkos::parallel_for(name, policy_fused_deidrj_y, *this);

    // z direction
    name[24] = '2';
    SnapAoSoATeamPolicy<DeviceType, team_size_compute_fused_deidrj, TagPairSNAPComputeFusedDeidrj<2> > policy_fused_deidrj_z(n_teams_div,team_size_compute_fused_deidrj,vector_length);
    policy_fused_deidrj_z = policy_fused_deidrj_z.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
    Kokkos::parallel_for(name, policy_fused_deidrj_z, *this);
  }

  // Helper routine that dispatches fully fused ComputeFusedDeidrj
  template <class TagPairSNAPComputeFusedDeidrjAll>
  void snap_dispatch_fused_deidrj_all(int n_teams_div, int scratch_size) {
    // make sure we're only passing in types we expect
    static_assert(std::is_same_v<TagPairSNAPComputeFusedDeidrjAllLarge, TagPairSNAPComputeFusedDeidrjAll> ||
      std::is_same_v<TagPairSNAPComputeFusedDeidrjAllSmall, TagPairSNAPComputeFusedDeidrjAll>);

    std::string name = { std::is_same_v<TagPairSNAPComputeFusedDeidrjAllLarge, TagPairSNAPComputeFusedDeidrjAll> ? "ComputeFusedDeidrjAllLarge" : "ComputeFusedDeidrjAllSmall" };

    // fully fused
    SnapAoSoATeamPolicy<DeviceType, team_size_compute_fused_deidrj, TagPairSNAPComputeFusedDeidrjAll> policy_fused_deidrj_all(n_teams_div,team_size_compute_fused_deidrj_all,vector_length);
    policy_fused_deidrj_all = policy_fused_deidrj_all.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
    Kokkos::parallel_for(name, policy_fused_deidrj_all, *this);
  }

  PairSNAPKokkos(class LAMMPS *);
  ~PairSNAPKokkos() override;

  void coeff(int, char**) override;
  void init_style() override;
  double init_one(int, int) override;
  void compute(int, int) override;
  double memory_usage() override;

  // CPU and GPU backend
  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const int& ii) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const int& ii, EV_FLOAT&) const;

  // GPU backend only
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeNeigh>::member_type& team) const;

  // GPU backend only
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeCayleyKlein, const int iatom_mod, const int jnbor, const int iatom_div) const;

  // CPU and GPU
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPPreUi, const int& iatom_mod, const int& j, const int& iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPPreUi, const int& iatom, const int& j) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPPreUi, const int& iatom) const;

  template<bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeUiSmall<chemsnap>,
    const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeUiSmall<chemsnap>>::member_type& team) const;

  template<bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeUiLarge<chemsnap>,
    const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeUiLarge<chemsnap>>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPTransformUi, const int& iatom_mod, const int& idxu, const int& iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPTransformUi, const int& iatom, const int& idxu) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPTransformUi, const int& iatom) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeZi<chemsnap>, const int& iatom_mod, const int& idxz, const int& iatom_div) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeZi<chemsnap>, const int& iatom, const int& idxz) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeZi<chemsnap>, const int& iatom) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBi<chemsnap>, const int& iatom_mod, const int& idxb, const int& iatom_div) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBi<chemsnap>, const int& iatom, const int& idxb) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBi<chemsnap>, const int& iatom) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBetaLinear, const int& iatom_mod, const int& idxb, const int& iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBetaLinear, const int& iatom, const int& idxb) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBetaLinear, const int& iatom) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBetaQuadratic, const int& iatom_mod, const int& idxb, const int& iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBetaQuadratic, const int& iatom, const int& idxb) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeBetaQuadratic, const int& iatom) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeYi<chemsnap>, const int& iatom_mod, const int& idxz, const int& iatom_div) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeYi<chemsnap>, const int& iatom, const int& idxz) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeYi<chemsnap>, const int& iatom) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeYiWithZlist<chemsnap>, const int& iatom_mod, const int& idxz, const int& iatom_div) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeYiWithZlist<chemsnap>, const int& iatom, const int& idxz) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeYiWithZlist<chemsnap>, const int& iatom) const;

  template<int dir> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeFusedDeidrjSmall<dir>,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeFusedDeidrjSmall<dir> >::member_type& team) const;

  template<int dir> KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeFusedDeidrjLarge<dir>,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeFusedDeidrjLarge<dir> >::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeFusedDeidrjAllSmall, const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeFusedDeidrjAllSmall>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeFusedDeidrjAllLarge, const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeFusedDeidrjAllLarge>::member_type& team) const;

  // CPU backend only
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeNeighCPU,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeNeighCPU>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeUiCPU, const int& iatom_mod, const int& idxu, const int& iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeUiCPU, const int& iatom, const int& jnbor) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeUiCPU, const int& iatom) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeDuidrjCPU, const int& iatom_mod, const int& jnbor, const int& iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeDuidrjCPU, const int& iatom, const int& jnbor) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeDuidrjCPU, const int& iatom) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeDeidrjCPU, const int& iatom_mod, const int& jnbor, const int& iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeDeidrjCPU, const int& iatom, const int& jnbor) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairSNAPComputeDeidrjCPU, const int& iatom) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const double &fx, const double &fy, const double &fz,
      const double &delx, const double &dely, const double &delz) const;

 protected:
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  SNAKokkos<DeviceType, real_type, vector_length> snaKK;

  int inum, max_neighs, batched_max_neighs, chunk_size, chunk_offset;
  int neighflag;

  int eflag,vflag;

  void allocate() override;

  Kokkos::View<real_type*, DeviceType> d_radelem;              // element radii
  Kokkos::View<real_type*, DeviceType> d_wjelem;               // elements weights
  typename SNAKokkos<DeviceType, real_type, vector_length>::t_sna_2d_lr d_coeffelem; // element bispectrum coefficients
  Kokkos::View<real_type*, DeviceType> d_sinnerelem;           // element inner cutoff midpoint
  Kokkos::View<real_type*, DeviceType> d_dinnerelem;           // element inner cutoff half-width
  Kokkos::View<T_INT*, DeviceType> d_map;                    // mapping from atom types to elements
  Kokkos::View<T_INT*, DeviceType> d_ninside;                // ninside for all atoms in list
  typename SNAKokkos<DeviceType, real_type, vector_length>::t_sna_2d d_beta;                // betas for all atoms in list

  typedef Kokkos::DualView<double**, DeviceType> tdual_fparams;
  tdual_fparams k_cutsq;
  typedef Kokkos::View<const double**, DeviceType,
      Kokkos::MemoryTraits<Kokkos::RandomAccess> > t_fparams_rnd;
  t_fparams_rnd rnd_cutsq;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;

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

  friend void pair_virial_fdotr_compute<PairSNAPKokkos>(PairSNAPKokkos*);

  // Utility routine which wraps computing per-team scratch size requirements for
  // ComputeNeigh, ComputeUi, and ComputeFusedDeidrj
  template <typename scratch_type>
  int scratch_size_helper(int values_per_team);

  // Make SNAKokkos a friend
  friend class SNAKokkos<DeviceType, real_type, vector_length>;
};


// These wrapper classes exist to make the pair style factory happy/avoid having
// to extend the pair style factory to support Pair classes w/an arbitrary number
// of extra template parameters

template <class DeviceType>
class PairSNAPKokkosDevice : public PairSNAPKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_DEVICE_VECLEN> {

 private:
  using Base = PairSNAPKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_DEVICE_VECLEN>;

 public:

  PairSNAPKokkosDevice(class LAMMPS *);

  void coeff(int, char**) override;
  void init_style() override;
  double init_one(int, int) override;
  void compute(int, int) override;
  double memory_usage() override;

};

#ifdef LMP_KOKKOS_GPU
template <class DeviceType>
class PairSNAPKokkosHost : public PairSNAPKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_HOST_VECLEN> {

 private:
  using Base = PairSNAPKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_HOST_VECLEN>;

 public:

  PairSNAPKokkosHost(class LAMMPS *);

  void coeff(int, char**);
  void init_style();
  double init_one(int, int);
  void compute(int, int);
  double memory_usage();

};
#endif

}

#endif
#endif

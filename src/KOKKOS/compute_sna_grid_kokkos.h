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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(sna/grid/kk,ComputeSNAGridKokkosDevice<LMPDeviceType>);
ComputeStyle(sna/grid/kk/device,ComputeSNAGridKokkosDevice<LMPDeviceType>);
#ifdef LMP_KOKKOS_GPU
ComputeStyle(sna/grid/kk/host,ComputeSNAGridKokkosHost<LMPHostType>);
#else
ComputeStyle(sna/grid/kk/host,ComputeSNAGridKokkosDevice<LMPHostType>);
#endif
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_SNA_GRID_KOKKOS_H
#define LMP_COMPUTE_SNA_GRID_KOKKOS_H

#include "compute_sna_grid.h"
#include "kokkos_type.h"
#include "sna_kokkos.h"

namespace LAMMPS_NS {

// Routines for both the CPU and GPU backend

// GPU backend only
struct TagCSNAGridComputeNeigh{};
struct TagCSNAGridComputeCayleyKlein{};
struct TagCSNAGridPreUi{};
struct TagCSNAGridComputeUiSmall{}; // more parallelism, more divergence
struct TagCSNAGridComputeUiLarge{}; // less parallelism, no divergence
struct TagCSNAGridTransformUi{}; // re-order ulisttot from SoA to AoSoA, zero ylist
template <bool chemsnap> struct TagCSNAGridComputeZi{};
template <bool chemsnap> struct TagCSNAGridComputeBi{};
struct TagCSNAGridLocalFill{}; // fill the gridlocal array

struct TagComputeSNAGridLoop{};
struct TagComputeSNAGrid3D{};

// CPU backend only
struct TagComputeSNAGridLoopCPU{};

//template<class DeviceType>
template<class DeviceType, typename real_type_, int vector_length_>
class ComputeSNAGridKokkos : public ComputeSNAGrid {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  static constexpr int vector_length = vector_length_;
  using real_type = real_type_;
  using complex = SNAComplex<real_type>;

  // Static team/tile sizes for device offload

#ifdef KOKKOS_ENABLE_HIP
  static constexpr int team_size_compute_neigh = 2;
  static constexpr int tile_size_compute_ck = 2;
  static constexpr int tile_size_pre_ui = 2;
  static constexpr int team_size_compute_ui = 2;
  static constexpr int tile_size_transform_ui = 2;
  static constexpr int tile_size_compute_zi = 2;
  static constexpr int min_blocks_compute_zi = 0; // no minimum bound
  static constexpr int tile_size_compute_bi = 2;
  static constexpr int tile_size_compute_yi = 2;
  static constexpr int min_blocks_compute_yi = 0; // no minimum bound
  static constexpr int team_size_compute_fused_deidrj = 2;
#else
  static constexpr int team_size_compute_neigh = 4;
  static constexpr int tile_size_compute_ck = 4;
  static constexpr int tile_size_pre_ui = 4;
  static constexpr int team_size_compute_ui = sizeof(real_type) == 4 ? 8 : 4;
  static constexpr int tile_size_transform_ui = 4;
  static constexpr int tile_size_compute_zi = 8;
  static constexpr int tile_size_compute_bi = 4;
  static constexpr int tile_size_compute_yi = 8;
  static constexpr int team_size_compute_fused_deidrj = sizeof(real_type) == 4 ? 4 : 2;

  // this empirically reduces perf fluctuations from compiler version to compiler version
  static constexpr int min_blocks_compute_zi = 4;
  static constexpr int min_blocks_compute_yi = 4;
#endif

  // Custom MDRangePolicy, Rank3, to reduce verbosity of kernel launches
  // This hides the Kokkos::IndexType<int> and Kokkos::Rank<3...>
  // and reduces the verbosity of the LaunchBound by hiding the explicit
  // multiplication by vector_length
  template <class Device, int num_tiles, class TagComputeSNA, int min_blocks = 0>
  using Snap3DRangePolicy = typename Kokkos::MDRangePolicy<Device, Kokkos::IndexType<int>, Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>, Kokkos::LaunchBounds<vector_length * num_tiles, min_blocks>, TagComputeSNA>;

  // MDRangePolicy for the 3D grid loop:
  template <class Device, class TagComputeSNA>
  using CSNAGrid3DPolicy = typename Kokkos::MDRangePolicy<Device, Kokkos::IndexType<int>, Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>;

  // Testing out team policies
  template <class Device, int num_teams,  class TagComputeSNA>
  using CSNAGridTeamPolicy = typename Kokkos::TeamPolicy<Device, Kokkos::LaunchBounds<vector_length * num_teams>, TagComputeSNA>;
  //using CSNAGridTeamPolicy = typename Kokkos::TeamPolicy<Device, Kokkos::IndexType<int>, Kokkos::IndexType<int>, Kokkos::IndexType<int>, TagComputeSNA>;
  //using team_member = typename team_policy::member_type;

  // Custom SnapAoSoATeamPolicy to reduce the verbosity of kernel launches
  // This hides the LaunchBounds abstraction by hiding the explicit
  // multiplication by vector length
  template <class Device, int num_teams, class TagComputeSNA>
  using SnapAoSoATeamPolicy = typename Kokkos::TeamPolicy<Device, Kokkos::LaunchBounds<vector_length * num_teams>, TagComputeSNA>;

  // Helper routine that returns a CPU or a GPU policy as appropriate
  template <class Device, int num_tiles, class TagComputeSNA, int min_blocks = 0>
  auto snap_get_policy(const int& chunk_size_div, const int& second_loop) {
    return Snap3DRangePolicy<Device, num_tiles, TagComputeSNA, min_blocks>({0, 0, 0},
                                                                 {vector_length, second_loop, chunk_size_div},
                                                                 {vector_length, num_tiles, 1});
  }

  ComputeSNAGridKokkos(class LAMMPS *, int, char **);
  ~ComputeSNAGridKokkos() override;

  void setup() override;
  void compute_array() override;

  // Utility functions for teams

  template<class TagStyle>
  void check_team_size_for(int, int&);

  template<class TagStyle>
  void check_team_size_reduce(int, int&);

  // operator function for example team policy
  //KOKKOS_INLINE_FUNCTION

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeSNAGridLoop, const int& ) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagComputeSNAGridLoopCPU, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagCSNAGridComputeNeigh>::member_type& team) const;

  // 3D case - used by parallel_for
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeSNAGrid3D, const int& iz, const int& iy, const int& ix) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridComputeCayleyKlein, const int iatom_mod, const int jnbor, const int iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridPreUi, const int& iatom_mod, const int& j, const int& iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridPreUi, const int& iatom, const int& j) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridPreUi, const int& iatom) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridComputeUiSmall,const typename Kokkos::TeamPolicy<DeviceType, TagCSNAGridComputeUiSmall>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridComputeUiLarge,const typename Kokkos::TeamPolicy<DeviceType, TagCSNAGridComputeUiLarge>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridTransformUi, const int& iatom_mod, const int& idxu, const int& iatom_div) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridTransformUi, const int& iatom, const int& idxu) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridTransformUi, const int& iatom) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridComputeZi<chemsnap>, const int& iatom_mod, const int& idxz, const int& iatom_div) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridComputeZi<chemsnap>, const int& iatom, const int& idxz) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridComputeZi<chemsnap>, const int& iatom) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridComputeBi<chemsnap>, const int& iatom_mod, const int& idxb, const int& iatom_div) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridComputeBi<chemsnap>, const int& iatom, const int& idxb) const;

  template <bool chemsnap> KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridComputeBi<chemsnap>, const int& iatom) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagCSNAGridLocalFill,const int& ii) const;

 protected:

  SNAKokkos<DeviceType, real_type, vector_length> snaKK;

  int max_neighs, chunk_size, chunk_offset;
  int host_flag;
  int ntotal;
  int total_range; // total number of loop iterations in grid
  int zlen; //= nzhi-nzlo+1;
  int ylen; //= nyhi-nylo+1;
  int xlen; //= nxhi-nxlo+1;

  double cutsq_tmp; // temporary cutsq until we get a view

  Kokkos::View<real_type*, DeviceType> d_radelem;              // element radii
  Kokkos::View<real_type*, DeviceType> d_wjelem;               // elements weights
  Kokkos::View<real_type**, Kokkos::LayoutRight, DeviceType> d_coeffelem;           // element bispectrum coefficients
  Kokkos::View<real_type*, DeviceType> d_sinnerelem;           // element inner cutoff midpoint
  Kokkos::View<real_type*, DeviceType> d_dinnerelem;           // element inner cutoff half-width
  Kokkos::View<T_INT*, DeviceType> d_ninside;                // ninside for all atoms in list
  Kokkos::View<T_INT*, DeviceType> d_map;                    // mapping from atom types to elements
  Kokkos::View<real_type*, DeviceType> d_test;              // test view

  typedef Kokkos::DualView<F_FLOAT**, DeviceType> tdual_fparams;
  tdual_fparams k_cutsq;
  typedef Kokkos::View<const F_FLOAT**, DeviceType,
      Kokkos::MemoryTraits<Kokkos::RandomAccess> > t_fparams_rnd;
  t_fparams_rnd rnd_cutsq;

  typename AT::t_x_array_randomread x;
  typename AT::t_int_1d_randomread type;
  DAT::tdual_float_2d k_grid;
  DAT::tdual_float_2d k_gridall;
  typename AT::t_float_2d d_grid;
  typename AT::t_float_2d d_gridall;

  DAT::tdual_float_4d k_gridlocal;
  typename AT::t_float_4d d_gridlocal;


  // Utility routine which wraps computing per-team scratch size requirements for
  // ComputeNeigh, ComputeUi, and ComputeFusedDeidrj
  template <typename scratch_type>
  int scratch_size_helper(int values_per_team);

  class DomainKokkos *domainKK;

  // triclinic vars
  double h0, h1, h2, h3, h4, h5;
  double lo0, lo1, lo2;

  // Make SNAKokkos a friend
  friend class SNAKokkos<DeviceType, real_type, vector_length>;
};

// These wrapper classes exist to make the compute style factory happy/avoid having
// to extend the compute  style factory to support Compute classes w/an arbitrary number
// of extra template parameters

template <class DeviceType>
class ComputeSNAGridKokkosDevice : public ComputeSNAGridKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_DEVICE_VECLEN> {

 private:
  using Base = ComputeSNAGridKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_DEVICE_VECLEN>;

 public:

  ComputeSNAGridKokkosDevice(class LAMMPS *, int, char **);

  void compute_array() override;

};

#ifdef LMP_KOKKOS_GPU
template <class DeviceType>
class ComputeSNAGridKokkosHost : public ComputeSNAGridKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_HOST_VECLEN> {

 private:
  using Base = ComputeSNAGridKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_HOST_VECLEN>;

 public:

  ComputeSNAGridKokkosHost(class LAMMPS *, int, char **);

  void compute_array() override;

};
#endif

}

#endif
#endif

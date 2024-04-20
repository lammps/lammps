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
FixStyle(wall/flow/kk,FixWallFlowKokkos<LMPDeviceType>);
FixStyle(wall/flow/kk/device,FixWallFlowKokkos<LMPDeviceType>);
FixStyle(wall/flow/kk/host,FixWallFlowKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_WALL_FLOW_KOKKOS_H
#define LMP_FIX_WALL_FLOW_KOKKOS_H

#include "fix_wall_flow.h"
#include "kokkos_type.h"
#include "kokkos_base.h"
#include "Kokkos_Random.hpp"

namespace LAMMPS_NS {

struct TagFixWallFlowInit{};
template<class MTag>
struct TagFixWallFlowEndOfStep{};
struct TagFixWallFlowPackExchange{};
struct TagFixWallFlowUnpackExchange{};

template<class DeviceType>
class FixWallFlowKokkos : public FixWallFlow, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  struct MassTag{};
  struct RMassTag{};
  FixWallFlowKokkos(class LAMMPS *, int, char **);
  ~FixWallFlowKokkos();

  void init() override;
  void end_of_step() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagFixWallFlowInit, const int&) const;

  template<class MTag>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallFlowEndOfStep<MTag>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallFlowPackExchange, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallFlowUnpackExchange, const int&) const;

  int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;

  void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              int /*nrecv1*/, int /*nextrarecv1*/,
                              ExecutionSpace space) override;
 protected:
  typename AT::t_x_array d_x;
  typename AT::t_v_array d_v;
  typename AT::t_int_1d d_type;
  typename AT::t_int_1d d_mask;

  typename AT::t_float_1d d_mass;
  typename AT::t_float_1d d_rmass;

  typedef typename AT::t_xfloat_1d d_walls_t;
  typedef Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool_t;
  typedef typename rand_pool_t::generator_type rand_type_t;

  typename AT::tdual_int_1d k_current_segment;
  typename AT::t_int_1d d_current_segment;
  typename HAT::t_int_1d h_current_segment;

  typename AT::t_int_1d d_sendlist;
  typename AT::t_xfloat_1d d_buf;
  typename AT::t_int_1d d_copylist;
  typename AT::t_int_1d d_indices;

  d_walls_t d_walls;

  rand_pool_t rand_pool;

  template<class MTag>
  KOKKOS_INLINE_FUNCTION
  void generate_velocity_kk(int atom_i) const;

  KOKKOS_INLINE_FUNCTION
  int compute_current_segment_kk(double pos) const;

  KOKKOS_INLINE_FUNCTION
  double get_mass(MassTag, int atom_i) const
  {
    return d_mass(d_type(atom_i));
  }

  KOKKOS_INLINE_FUNCTION
  double get_mass(RMassTag, int atom_i) const
  {
    return d_rmass(atom_i);
  }
};

}

#endif
#endif


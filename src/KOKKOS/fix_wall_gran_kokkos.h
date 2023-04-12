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
FixStyle(wall/gran/kk,FixWallGranKokkos<LMPDeviceType>)
FixStyle(wall/gran/kk/device,FixWallGranKokkos<LMPDeviceType>)
FixStyle(wall/gran/kk/host,FixWallGranKokkos<LMPHostType>)
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_WALL_GRAN_KOKKOS_H
#define LMP_FIX_WALL_GRAN_KOKKOS_H

#include "fix_wall_gran_old.h"
#include "kokkos_type.h"
#include "kokkos_base.h"

namespace LAMMPS_NS {

template<int WallStyle>
struct TagFixWallGranHookeHistory{};

struct TagFixWallGranPackExchange{};
struct TagFixWallGranUnpackExchange{};

template<class DeviceType>
class FixWallGranKokkos : public FixWallGranOld, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixWallGranKokkos(class LAMMPS *, int, char **);
  ~FixWallGranKokkos() override;
  void init() override;
  void post_force(int) override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  template <int WallStyle>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallGranHookeHistory<WallStyle>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallGranPackExchange, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallGranUnpackExchange, const int&) const;

  int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
			   DAT::tdual_int_1d k_sendlist,
			   DAT::tdual_int_1d k_copylist,
			   ExecutionSpace space) override;

  void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              ExecutionSpace space) override;

 private:
  X_FLOAT wlo;
  X_FLOAT whi;
  V_FLOAT vwall[3];

  typename AT::t_x_array x;
  typename AT::t_v_array v;
  typename AT::t_v_array d_omega;
  typename AT::t_f_array f;
  typename AT::t_f_array torque;
  typename AT::t_int_1d mask;
  typename AT::t_float_1d rmass;
  typename AT::t_float_1d d_radius;
  typename AT::tdual_float_2d k_history_one;
  typename AT::t_float_2d d_history_one;

  typename AT::t_int_1d d_sendlist;
  typename AT::t_xfloat_1d d_buf;
  typename AT::t_int_1d d_copylist;
  typename AT::t_int_1d d_indices;
};
}

#endif
#endif

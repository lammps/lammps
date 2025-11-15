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
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  template <int WallStyle>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallGranHookeHistory<WallStyle>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallGranPackExchange, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallGranUnpackExchange, const int&) const;

  int pack_exchange_kokkos(const int &nsend,DAT::tdual_double_2d_lr &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;

  void unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              int nrecv1,int nrecv1extra,
                              ExecutionSpace space) override;

 private:
  KK_FLOAT wlo;
  KK_FLOAT whi;
  KK_FLOAT vwall[3];

  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkfloat_1d_3 d_omega;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d mask;
  typename AT::t_kkfloat_1d rmass;
  typename AT::t_kkfloat_1d d_radius;
  DAT::ttransform_kkfloat_2d k_history_one;
  typename AT::t_kkfloat_2d d_history_one;

  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d d_buf;
  typename AT::t_int_1d d_copylist;
  typename AT::t_int_1d d_indices;
};

}

#endif
#endif

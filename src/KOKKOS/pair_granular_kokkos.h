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
PairStyle(granular/kk,PairGranularKokkos<LMPDeviceType>);
PairStyle(granular/kk/device,PairGranularKokkos<LMPDeviceType>);
PairStyle(granular/kk/host,PairGranularKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_GRANULAR_KOKKOS_H
#define LMP_PAIR_GRANULAR_KOKKOS_H

#include "kokkos_type.h"
#include "pair_kokkos.h"

#include "pair_granular.h"

// the shared sub-model kernels; kokkos_type.h above makes them device callable
#include "gran_sub_mod_damping_kernel.h"
#include "gran_sub_mod_normal_kernel.h"
#include "gran_sub_mod_rolling_kernel.h"
#include "gran_sub_mod_tangential_kernel.h"
#include "gran_sub_mod_twisting_kernel.h"

namespace LAMMPS_NS {

template <class DeviceType> class FixNeighHistoryKokkos;

// largest per-contact history of any supported combination:
// tangential (4) + rolling (3) + twisting (3)
static constexpr int GRAN_KK_MAX_HISTORY = 10;

// per sub-model history held in registers inside the kernel.  every index
// into these arrays must be a compile-time constant, or the whole array
// goes to local memory, which costs ~1.4x in kernel time on a GP100.
static constexpr int GRAN_KK_TANGENTIAL_HISTORY = 4;
static constexpr int GRAN_KK_ROLLING_HISTORY = 3;
static constexpr int GRAN_KK_TWISTING_HISTORY = 3;

// EXTRA selects whether the rolling and twisting models are compiled in.
// Leaving them out of the common case keeps the register footprint down.

template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int EXTRA>
struct TagPairGranularCompute {};

// everything the kernel needs to know about one pair_coeff model

template <class T>
struct GranularModelKK {
  Granular_NS::GranKernel::GranNormalParams<T> normal;
  Granular_NS::GranKernel::GranDampingParams<T> damping;
  Granular_NS::GranKernel::GranTangentialParams<T> tangential;
  Granular_NS::GranKernel::GranRollingParams<T> rolling;
  Granular_NS::GranKernel::GranTwistingParams<T> twisting;
  int tangential_index;
  int rolling_index;
  int twisting_index;
  int tangential_size;
  int rolling_size;
  int twisting_size;
  int limit_damping;
  int contact_radius_flag;
};

template <class DeviceType>
class PairGranularKokkos : public PairGranular {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairGranularKokkos(class LAMMPS *);
  ~PairGranularKokkos() override;
  void compute(int, int) override;
  void init_style() override;
  void setup() override;
  std::string history_fix_command() override;

  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int EXTRA>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranularCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,EXTRA>, const int, EV_FLOAT &) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int EXTRA>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranularCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,EXTRA>, const int) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, int i, int j,
                    KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                    KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const;

 protected:
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_randomread v;
  typename AT::t_kkfloat_1d_3_randomread omega;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread radius;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  typename AT::t_int_2d d_firsttouch;
  typename AT::t_kkfloat_2d d_firsthistory;

  Kokkos::View<GranularModelKK<KK_FLOAT>*,DeviceType> d_models;
  Kokkos::View<int**,Kokkos::LayoutRight,DeviceType> d_type_model;

  FixNeighHistoryKokkos<DeviceType> *fix_historyKK;

  int neighflag;
  int newton_pair;
  int nlocal,nall,eflag,vflag;
  int history_update;
  int size_history_kk;
  int use_history_kk;
  int extra_models;         // rolling and/or twisting present
  KK_FLOAT dt_kk;
  KK_FLOAT special_lj[4];

  void create_history_fix();
  void pack_models();

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const { return j >> SBBITS & 3; }

  friend void pair_virial_fdotr_compute<PairGranularKokkos>(PairGranularKokkos*);
};

}    // namespace LAMMPS_NS

#endif
#endif

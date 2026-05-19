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
FixStyle(gjf/kk,FixGJFKokkos<LMPDeviceType>);
FixStyle(gjf/kk/device,FixGJFKokkos<LMPDeviceType>);
FixStyle(gjf/kk/host,FixGJFKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_GJF_KOKKOS_H
#define LMP_FIX_GJF_KOKKOS_H

#include "fix_gjf.h"
#include "kokkos_base.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap_kokkos.h"

namespace LAMMPS_NS {

struct TagFixGJFInitialIntegrate{};
struct TagFixGJFFinalIntegrate{};
struct TagFixGJFEndOfStep{};
struct TagFixGJFUnpackExchange{};

template<class DeviceType>
class FixGJFKokkos : public FixGJF, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixGJFKokkos(class LAMMPS *, int, char **);
  ~FixGJFKokkos() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void end_of_step() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixGJFInitialIntegrate, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixGJFFinalIntegrate, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixGJFEndOfStep, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void pack_exchange_item(const int &, int &, const bool &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixGJFUnpackExchange, const int &) const;

  int pack_exchange_kokkos(const int &nsend, DAT::tdual_double_2d_lr &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;

  void unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf,
                              DAT::tdual_int_1d &indices, int nrecv,
                              int nrecv1, int nrecv1extra,
                              ExecutionSpace space) override;

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

 protected:
  int nrecv1, nextrarecv1;

  // per-atom half-step velocity storage (replaces base class lv[][3])
  DAT::ttransform_kkfloat_1d_3_lr k_lv;
  typename AT::t_kkfloat_1d_3_lr d_lv;

  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkfloat_1d_3_lr v;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;

  int nsend;
  typename AT::t_int_2d d_sendlist;
  typename AT::t_double_1d_um d_buf;
  typename AT::t_int_1d d_exchange_sendlist;
  typename AT::t_int_1d d_copylist;
  typename AT::t_int_1d d_indices;

  typename AT::t_int_scalar d_count;
  HAT::t_int_scalar h_count;

  // cached integration constants (set in initial_integrate / final_integrate)
  KK_FLOAT l_gjfc2, l_c1sqrt, l_c3sqrt, l_csq;
  KK_FLOAT l_dtf, l_dt, l_tsqrt, l_t_period;
  KK_FLOAT l_boltz, l_mvv2e, l_ftm2v;
  int l_rmass_flag;

#ifndef LMP_KOKKOS_DEBUG_RNG
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif
};

template<class DeviceType>
struct FixGJFKokkosPackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;
  FixGJFKokkos<DeviceType> c;
  FixGJFKokkosPackExchangeFunctor(FixGJFKokkos<DeviceType>* c_ptr) : c(*c_ptr) {}
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i, int &offset, const bool &final) const {
    c.pack_exchange_item(i, offset, final);
  }
};

}    // namespace LAMMPS_NS

#endif
#endif

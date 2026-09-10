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
FixStyle(baoab/kk,FixBAOABKokkos<LMPDeviceType>);
FixStyle(baoab/kk/device,FixBAOABKokkos<LMPDeviceType>);
FixStyle(baoab/kk/host,FixBAOABKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_BAOAB_KOKKOS_H
#define LMP_FIX_BAOAB_KOKKOS_H

#include "fix_baoab.h"

#include "kokkos_type.h"
#include "rand_pool_wrap_kokkos.h"

#include <Kokkos_Random.hpp>

namespace LAMMPS_NS {

// RMASS selects per-atom vs per-type mass; ZERO enables the accumulation
// needed by the zero-total-random-momentum correction

template<int RMASS, int ZERO> struct TagFixBAOABInitial{};
template<int RMASS> struct TagFixBAOABFinal{};
struct TagFixBAOABZeroMomentum{};

template<class DeviceType>
class FixBAOABKokkos : public FixBAOAB {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];
  const int value_count = 5;

  FixBAOABKokkos(class LAMMPS *, int, char **);
  ~FixBAOABKokkos() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;

  template<int RMASS, int ZERO>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixBAOABInitial<RMASS,ZERO>, const int&, value_type) const;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixBAOABFinal<RMASS>, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixBAOABZeroMomentum, const int&) const;

 private:
  typename AT::t_kkfloat_1d_3_lr d_x;
  typename AT::t_kkfloat_1d_3 d_v;
  typename AT::t_kkacc_1d_3_randomread d_f;
  typename AT::t_kkfloat_1d_randomread d_rmass;
  typename AT::t_kkfloat_1d_randomread d_mass;
  typename AT::t_int_1d_randomread d_type;
  typename AT::t_int_1d_randomread d_mask;

  KK_FLOAT l_dtf, l_dtby2, l_c1, l_kT, l_one_minus_c1sq, l_mvv2e;
  KK_FLOAT l_vcorr[3];

#ifndef LMP_KOKKOS_DEBUG_RNG
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif
};

}

#endif
#endif

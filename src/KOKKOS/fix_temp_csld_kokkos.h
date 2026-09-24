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
FixStyle(temp/csld/kk,FixTempCSLDKokkos<LMPDeviceType>);
FixStyle(temp/csld/kk/device,FixTempCSLDKokkos<LMPDeviceType>);
FixStyle(temp/csld/kk/host,FixTempCSLDKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_TEMP_CSLD_KOKKOS_H
#define LMP_FIX_TEMP_CSLD_KOKKOS_H

#include "fix_temp_csld.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap_kokkos.h"

namespace LAMMPS_NS {

struct TagFixTempCSLDInitial{};
struct TagFixTempCSLDFinal{};

template<class DeviceType>
class FixTempCSLDKokkos : public FixTempCSLD {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixTempCSLDKokkos(class LAMMPS *, int, char **);
  ~FixTempCSLDKokkos() override;
  void init() override;
  void end_of_step() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixTempCSLDInitial, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixTempCSLDFinal, const int&) const;

 protected:
  double temp_scalar();

  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkfloat_1d_3 d_vhold;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;

  int l_rmass_flag;
  KK_FLOAT l_c1, l_c2;

#ifndef LMP_KOKKOS_DEBUG_RNG
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif
};

}    // namespace LAMMPS_NS

#endif
#endif

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
FixStyle(brownian/kk,FixBrownianKokkos<LMPDeviceType>);
FixStyle(brownian/kk/device,FixBrownianKokkos<LMPDeviceType>);
FixStyle(brownian/kk/host,FixBrownianKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_BROWNIAN_KOKKOS_H
#define LMP_FIX_BROWNIAN_KOKKOS_H

#include "fix_brownian.h"

#include "kokkos_type.h"
#include "rand_pool_wrap_kokkos.h"

#include <Kokkos_Random.hpp>

namespace LAMMPS_NS {

// the template parameters mirror the CPU style's initial_integrate_templated()

template<int Tp_UNIFORM, int Tp_GAUSS, int Tp_2D>
struct TagFixBrownian{};

template<class DeviceType>
class FixBrownianKokkos : public FixBrownian {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixBrownianKokkos(class LAMMPS *, int, char **);
  ~FixBrownianKokkos() override;
  void init() override;
  void initial_integrate(int) override;

  template<int Tp_UNIFORM, int Tp_GAUSS, int Tp_2D>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixBrownian<Tp_UNIFORM,Tp_GAUSS,Tp_2D>, const int &) const;

 private:
  template<int Tp_UNIFORM, int Tp_GAUSS, int Tp_2D> void initial_integrate_kokkos();

  typename AT::t_kkfloat_1d_3_lr d_x;
  typename AT::t_kkfloat_1d_3 d_v;
  typename AT::t_kkacc_1d_3_randomread d_f;
  typename AT::t_int_1d_randomread d_mask;

  KK_FLOAT l_dt, l_g1, l_g2;

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

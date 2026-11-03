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
FixStyle(flow/gauss/kk,FixFlowGaussKokkos<LMPDeviceType>);
FixStyle(flow/gauss/kk/device,FixFlowGaussKokkos<LMPDeviceType>);
FixStyle(flow/gauss/kk/host,FixFlowGaussKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_FLOW_GAUSS_KOKKOS_H
#define LMP_FIX_FLOW_GAUSS_KOKKOS_H

#include "fix_flow_gauss.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

// three-component force sum; a class-level "typedef double value_type[]"
// would force the added-work reduction below into the array form too

struct s_KK_double3 {
  double d0, d1, d2;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  s_KK_double3() {
    d0 = d1 = d2 = 0.0;
  }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  s_KK_double3& operator+=(const s_KK_double3 &rhs) {
    d0 += rhs.d0;
    d1 += rhs.d1;
    d2 += rhs.d2;
    return *this;
  }
};

struct TagFixFlowGaussReduce{};
struct TagFixFlowGaussApply{};
struct TagFixFlowGaussApplyWork{};

template<class DeviceType>
class FixFlowGaussKokkos : public FixFlowGauss {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixFlowGaussKokkos(class LAMMPS *, int, char **);
  void init() override;
  void setup(int) override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixFlowGaussReduce, const int &, s_KK_double3 &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixFlowGaussApply, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixFlowGaussApplyWork, const int &, double &) const;

 private:
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_kkfloat_1d_3_randomread d_v;
  typename AT::t_int_1d_randomread d_mask;
  typename AT::t_int_1d_randomread d_type;
  typename AT::t_kkfloat_1d_randomread d_rmass;
  typename AT::t_kkfloat_1d_randomread d_mass;

  KK_FLOAT m_a_app[3];

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void applied_force(const int &i, KK_FLOAT *f_app) const
  {
    const KK_FLOAT massone = (d_rmass.data() ? d_rmass(i) : d_mass(d_type(i)));
    f_app[0] = m_a_app[0]*massone;
    f_app[1] = m_a_app[1]*massone;
    f_app[2] = m_a_app[2]*massone;
  }
};

}

#endif
#endif

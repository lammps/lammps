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
FixStyle(setforce/kk,FixSetForceKokkos<LMPDeviceType>);
FixStyle(setforce/kk/device,FixSetForceKokkos<LMPDeviceType>);
FixStyle(setforce/kk/host,FixSetForceKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_SET_FORCE_KOKKOS_H
#define LMP_FIX_SET_FORCE_KOKKOS_H

#include "fix_setforce.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct s_double_3 {
  double d0, d1, d2;
  KOKKOS_INLINE_FUNCTION
  s_double_3() {
    d0 = d1 = d2 = 0.0;
  }
  KOKKOS_INLINE_FUNCTION
  s_double_3& operator+=(const s_double_3 &rhs) {
    d0 += rhs.d0;
    d1 += rhs.d1;
    d2 += rhs.d2;
    return *this;
  }
};
typedef s_double_3 double_3;

struct TagFixSetForceConstant{};

struct TagFixSetForceNonConstant{};

template<class DeviceType>
class FixSetForceKokkos : public FixSetForce {
 public:
  typedef DeviceType device_type;
  typedef double_3 value_type;
  typedef ArrayTypes<DeviceType> AT;

  FixSetForceKokkos(class LAMMPS *, int, char **);
  ~FixSetForceKokkos() override;
  void init() override;
  void post_force(int) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSetForceConstant, const int&, double_3&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSetForceNonConstant, const int&, double_3&) const;

 private:
  DAT::tdual_ffloat_2d k_sforce;
  typename AT::t_ffloat_2d_randomread d_sforce;
  typename AT::t_int_1d d_match;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread mask;
};

}

#endif
#endif


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
FixStyle(efield/kk,FixEfieldKokkos<LMPDeviceType>);
FixStyle(efield/kk/device,FixEfieldKokkos<LMPDeviceType>);
FixStyle(efield/kk/host,FixEfieldKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_EFIELD_KOKKOS_H
#define LMP_FIX_EFIELD_KOKKOS_H

#include "fix_efield.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct e_double_4 {
  double d0, d1, d2, d3;
  KOKKOS_INLINE_FUNCTION
  e_double_4() {
    d0 = d1 = d2 = d3 = 0.0;
  }
  KOKKOS_INLINE_FUNCTION
  e_double_4& operator+=(const e_double_4 &rhs) {
    d0 += rhs.d0;
    d1 += rhs.d1;
    d2 += rhs.d2;
    d3 += rhs.d3;
    return *this;
  }
};
typedef e_double_4 double_4;

struct TagFixEfieldConstant{};

struct TagFixEfieldNonConstant{};

template<class DeviceType>
class FixEfieldKokkos : public FixEfield {
 public:
  typedef DeviceType device_type;
  typedef double_4 value_type;
  typedef ArrayTypes<DeviceType> AT;

  FixEfieldKokkos(class LAMMPS *, int, char **);
  ~FixEfieldKokkos() override;
  void init() override;
  void post_force(int) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEfieldConstant, const int&, double_4&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEfieldNonConstant, const int&, double_4&) const;

 private:
  DAT::tdual_ffloat_2d k_efield;
  typename AT::t_ffloat_2d_randomread d_efield;
  typename AT::t_int_1d d_match;

  typename AT::t_x_array_randomread x;
  typename AT::t_float_1d_randomread q;
  typename AT::t_f_array f;
  typename AT::t_imageint_1d_randomread image;
  typename AT::t_int_1d_randomread mask;
};

}

#endif
#endif


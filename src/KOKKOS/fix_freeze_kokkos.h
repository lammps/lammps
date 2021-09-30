/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(freeze/kk,FixFreezeKokkos<LMPDeviceType>);
FixStyle(freeze/kk/device,FixFreezeKokkos<LMPDeviceType>);
FixStyle(freeze/kk/host,FixFreezeKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_FREEZE_KOKKOS_H
#define LMP_FIX_FREEZE_KOKKOS_H

#include "fix_freeze.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixFreezeKokkos : public FixFreeze {
 public:
  typedef DeviceType device_type;
  struct OriginalForce {
    double values[3];

    KOKKOS_INLINE_FUNCTION
    OriginalForce() {
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
    }

    KOKKOS_INLINE_FUNCTION
    OriginalForce &operator+=(const OriginalForce &rhs) {
      values[0] += rhs.values[0];
      values[1] += rhs.values[1];
      values[2] += rhs.values[2];
      return *this;
    }

    KOKKOS_INLINE_FUNCTION
    void operator+=(const volatile OriginalForce &rhs) volatile {
      values[0] += rhs.values[0];
      values[1] += rhs.values[1];
      values[2] += rhs.values[2];
    }
  };

  FixFreezeKokkos(class LAMMPS *, int, char **);
  void post_force(int);

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, OriginalForce &original) const;

 private:
  typename ArrayTypes<DeviceType>::t_f_array f;
  typename ArrayTypes<DeviceType>::t_f_array torque;
  typename ArrayTypes<DeviceType>::t_int_1d mask;
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_FREEZE_KOKKOS_H
#endif // FIX_CLASS

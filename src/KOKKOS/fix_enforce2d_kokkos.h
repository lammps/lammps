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
FixStyle(enforce2d/kk,FixEnforce2DKokkos<LMPDeviceType>);
FixStyle(enforce2d/kk/device,FixEnforce2DKokkos<LMPDeviceType>);
FixStyle(enforce2d/kk/host,FixEnforce2DKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_ENFORCE2D_KOKKOS_H
#define LMP_FIX_ENFORCE2D_KOKKOS_H

#include "fix_enforce2d.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixEnforce2DKokkos : public FixEnforce2D {
 public:
  typedef ArrayTypes<DeviceType> AT;

  FixEnforce2DKokkos(class LAMMPS *, int, char **);
  // ~FixEnforce2DKokkos() {}
  void setup(int) override;
  void post_force(int) override;

  template <int OMEGA_FLAG, int ANGMOM_FLAG, int TORQUE_FLAG>
  KOKKOS_INLINE_FUNCTION
  void post_force_item(const int i) const;

  // void min_setup(int);       Kokkos does not support minimization (yet)
  // void min_post_force(int);  Kokkos does not support minimization (yet)
  // void post_force_respa(int, int, int);  No RRESPA support yet.

 private:
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkacc_1d_3 f;

  typename AT::t_kkfloat_1d_3 omega;
  typename AT::t_kkfloat_1d_3 angmom;
  typename AT::t_kkacc_1d_3 torque;

  typename AT::t_int_1d mask;
};


template <class DeviceType, int omega_flag, int angmom_flag, int torque_flag>
struct FixEnforce2DKokkosPostForceFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  FixEnforce2DKokkos<DeviceType> c;

  FixEnforce2DKokkosPostForceFunctor(FixEnforce2DKokkos<DeviceType>* c_ptr):
    c(*c_ptr) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    c.template post_force_item <omega_flag, angmom_flag, torque_flag>(i);
  }
};


}

#endif
#endif


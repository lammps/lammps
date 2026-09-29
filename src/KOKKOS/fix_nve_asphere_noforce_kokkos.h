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
FixStyle(nve/asphere/noforce/kk,FixNVEAsphereNoforceKokkos<LMPDeviceType>);
FixStyle(nve/asphere/noforce/kk/device,FixNVEAsphereNoforceKokkos<LMPDeviceType>);
FixStyle(nve/asphere/noforce/kk/host,FixNVEAsphereNoforceKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_NVE_ASPHERE_NOFORCE_KOKKOS_H
#define LMP_FIX_NVE_ASPHERE_NOFORCE_KOKKOS_H

#include "fix_nve_asphere_noforce.h"

#include "atom_vec_ellipsoid_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixNVEAsphereNoforceKokkos : public FixNVEAsphereNoforce {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixNVEAsphereNoforceKokkos(class LAMMPS *, int, char **);

  void cleanup_copy();
  void init() override;
  void initial_integrate(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void initial_integrate_item(const int i) const;

 private:
  class AtomVecEllipsoidKokkos *avecEllipKK;
  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d bonus;
  typename AT::t_int_1d ellipsoid;
  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkfloat_1d_3 angmom;
  typename AT::t_kkfloat_1d rmass;
  typename AT::t_int_1d mask;
};

template<class DeviceType>
struct FixNVEAsphereNoforceKokkosInitialIntegrateFunctor {
  typedef DeviceType device_type;
  FixNVEAsphereNoforceKokkos<DeviceType> c;
  FixNVEAsphereNoforceKokkosInitialIntegrateFunctor(FixNVEAsphereNoforceKokkos<DeviceType> *c_ptr):
    c(*c_ptr) { c.cleanup_copy(); }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    c.initial_integrate_item(i);
  }
};

}

#endif
#endif

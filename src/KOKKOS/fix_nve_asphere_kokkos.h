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
FixStyle(nve/asphere/kk,FixNVEAsphereKokkos<LMPDeviceType>);
FixStyle(nve/asphere/kk/device,FixNVEAsphereKokkos<LMPDeviceType>);
FixStyle(nve/asphere/kk/host,FixNVEAsphereKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_NVE_ASPHERE_KOKKOS_H
#define LMP_FIX_NVE_ASPHERE_KOKKOS_H

#include "fix_nve_asphere.h"
#include "kokkos_type.h"

#include "atom_vec_ellipsoid_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixNVEAsphereKokkos : public FixNVEAsphere {
  public:
    FixNVEAsphereKokkos(class LAMMPS *, int, char **);

    void cleanup_copy();
    void init() override;
    void initial_integrate(int) override;
    void final_integrate() override;
    void fused_integrate(int) override;

    KOKKOS_INLINE_FUNCTION
    void initial_integrate_item(const int i) const;
    KOKKOS_INLINE_FUNCTION
    void final_integrate_item(const int i) const;
    KOKKOS_INLINE_FUNCTION
    void fused_integrate_item(int) const;

  private:
    class AtomVecEllipsoidKokkos *avecEllipKK;
    typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d bonus;
    typename ArrayTypes<DeviceType>::t_int_1d ellipsoid;
    typename ArrayTypes<DeviceType>::t_kkfloat_1d_3_lr x;
    typename ArrayTypes<DeviceType>::t_kkfloat_1d_3 v;
    typename ArrayTypes<DeviceType>::t_kkacc_1d_3 f;
    typename ArrayTypes<DeviceType>::t_kkfloat_1d_3 angmom;
    typename ArrayTypes<DeviceType>::t_kkacc_1d_3 torque;
    typename ArrayTypes<DeviceType>::t_kkfloat_1d rmass;
    typename ArrayTypes<DeviceType>::t_int_1d mask;
};

template <class DeviceType>
struct FixNVEAsphereKokkosInitialIntegrateFunctor {
  typedef DeviceType device_type;
  FixNVEAsphereKokkos<DeviceType> c;
  FixNVEAsphereKokkosInitialIntegrateFunctor(FixNVEAsphereKokkos<DeviceType> *c_ptr): c(*c_ptr) { c.cleanup_copy(); }
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    c.initial_integrate_item(i);
  }
};

template <class DeviceType>
struct FixNVEAsphereKokkosFinalIntegrateFunctor {
  typedef DeviceType device_type;
  FixNVEAsphereKokkos<DeviceType> c;
  FixNVEAsphereKokkosFinalIntegrateFunctor(FixNVEAsphereKokkos<DeviceType> *c_ptr): c(*c_ptr) { c.cleanup_copy(); }
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    c.final_integrate_item(i);
  }
};

template <class DeviceType>
struct FixNVEAsphereKokkosFusedIntegrateFunctor {
  typedef DeviceType device_type;
  FixNVEAsphereKokkos<DeviceType> c;
  FixNVEAsphereKokkosFusedIntegrateFunctor(FixNVEAsphereKokkos<DeviceType> *c_ptr): c(*c_ptr) { c.cleanup_copy(); }
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    c.fused_integrate_item(i);
  }
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_NVE_ASPHERE_KOKKOS_H
#endif // FIX_CLASS

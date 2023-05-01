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
FixStyle(nve/kk,FixNVEKokkos<LMPDeviceType>);
FixStyle(nve/kk/device,FixNVEKokkos<LMPDeviceType>);
FixStyle(nve/kk/host,FixNVEKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_NVE_KOKKOS_H
#define LMP_FIX_NVE_KOKKOS_H

#include "fix_nve.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixNVEKokkos;

template <class DeviceType, int RMass>
struct FixNVEKokkosInitialIntegrateFunctor;

template <class DeviceType, int RMass>
struct FixNVEKokkosFinalIntegrateFunctor;

template<class DeviceType>
class FixNVEKokkos : public FixNVE {
 public:
  FixNVEKokkos(class LAMMPS *, int, char **);

  void cleanup_copy();
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void fused_integrate(int) override;

  KOKKOS_INLINE_FUNCTION
  void initial_integrate_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void initial_integrate_rmass_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void final_integrate_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void final_integrate_rmass_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void fused_integrate_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void fused_integrate_rmass_item(int) const;

 private:


  typename ArrayTypes<DeviceType>::t_x_array x;
  typename ArrayTypes<DeviceType>::t_v_array v;
  typename ArrayTypes<DeviceType>::t_f_array_const f;
  typename ArrayTypes<DeviceType>::t_float_1d rmass;
  typename ArrayTypes<DeviceType>::t_float_1d mass;
  typename ArrayTypes<DeviceType>::t_int_1d type;
  typename ArrayTypes<DeviceType>::t_int_1d mask;
};

template <class DeviceType, int RMass>
struct FixNVEKokkosInitialIntegrateFunctor  {
  typedef DeviceType  device_type ;
  FixNVEKokkos<DeviceType> c;

  FixNVEKokkosInitialIntegrateFunctor(FixNVEKokkos<DeviceType>* c_ptr):
  c(*c_ptr) {c.cleanup_copy();};
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (RMass) c.initial_integrate_rmass_item(i);
    else c.initial_integrate_item(i);
  }
};

template <class DeviceType, int RMass>
struct FixNVEKokkosFinalIntegrateFunctor  {
  typedef DeviceType  device_type ;
  FixNVEKokkos<DeviceType> c;

  FixNVEKokkosFinalIntegrateFunctor(FixNVEKokkos<DeviceType>* c_ptr):
  c(*c_ptr) {c.cleanup_copy();};
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (RMass) c.final_integrate_rmass_item(i);
    else c.final_integrate_item(i);
  }
};

template <class DeviceType, int RMass>
struct FixNVEKokkosFusedIntegrateFunctor  {
  typedef DeviceType  device_type ;
  FixNVEKokkos<DeviceType> c;

  FixNVEKokkosFusedIntegrateFunctor(FixNVEKokkos<DeviceType>* c_ptr):
  c(*c_ptr) {c.cleanup_copy();};
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (RMass)
      c.fused_integrate_rmass_item(i);
    else
      c.fused_integrate_item(i);
  }
};

}

#endif
#endif


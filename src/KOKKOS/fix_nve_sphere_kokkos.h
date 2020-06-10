/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nve/sphere/kk,FixNVESphereKokkos<Device>)
FixStyle(nve/sphere/kk/device,FixNVESphereKokkos<Device>)
FixStyle(nve/sphere/kk/host,FixNVESphereKokkos<Host>)

#else

#ifndef LMP_FIX_NVE_SPHERE_KOKKOS_H
#define LMP_FIX_NVE_SPHERE_KOKKOS_H

#include "fix_nve_sphere.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class FixNVESphereKokkos : public FixNVESphere {
  public:
    typedef typename GetDeviceType<Space>::value DeviceType;
    typedef ArrayTypes<Space> AT;

    FixNVESphereKokkos(class LAMMPS *, int, char **);
    virtual ~FixNVESphereKokkos() {}
    void cleanup_copy();
    void init();
    void initial_integrate(int);
    void final_integrate();

    KOKKOS_INLINE_FUNCTION
    void initial_integrate_item(const int i) const;
    KOKKOS_INLINE_FUNCTION
    void final_integrate_item(const int i) const;

  private:
    typename AT::t_float_1d_3_lr x;
    typename AT::t_float_1d_3 v;
    typename AT::t_float_1d_3 omega;
    typename AT::t_float_1d_3 f;
    typename AT::t_float_1d_3 torque;
    typename AT::t_float_1d rmass;
    typename AT::t_float_1d radius;
    typename AT::t_int_1d mask;
};

template <ExecutionSpace Space>
struct FixNVESphereKokkosInitialIntegrateFunctor {
  FixNVESphereKokkos<Space> c;
  FixNVESphereKokkosInitialIntegrateFunctor(FixNVESphereKokkos<Space> *c_ptr): c(*c_ptr) { c.cleanup_copy(); }
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    c.initial_integrate_item(i);
  }
};

template <ExecutionSpace Space>
struct FixNVESphereKokkosFinalIntegrateFunctor {
  FixNVESphereKokkos<Space> c;
  FixNVESphereKokkosFinalIntegrateFunctor(FixNVESphereKokkos<Space> *c_ptr): c(*c_ptr) { c.cleanup_copy(); }
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    c.final_integrate_item(i);
  }
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_NVE_SPHERE_KOKKOS_H
#endif // FIX_CLASS

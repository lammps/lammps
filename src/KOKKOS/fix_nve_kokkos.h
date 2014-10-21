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

FixStyle(nve/kk,FixNVEKokkos<LMPDeviceType>)
FixStyle(nve/kk/device,FixNVEKokkos<LMPDeviceType>)
FixStyle(nve/kk/host,FixNVEKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_NVE_KOKKOS_H
#define LMP_FIX_NVE_KOKKOS_H

#include "fix_nve.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixNVEKokkos;

template <class DeviceType, int RMass>
class FixNVEKokkosInitialIntegrateFunctor;
template <class DeviceType, int RMass>
class FixNVEKokkosFinalIntegrateFunctor;

template<class DeviceType>
class FixNVEKokkos : public FixNVE {
 public:
  FixNVEKokkos(class LAMMPS *, int, char **);
  ~FixNVEKokkos() {}
  void cleanup_copy();
  void init();
  void initial_integrate(int);
  void final_integrate();

  KOKKOS_INLINE_FUNCTION
  void initial_integrate_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void initial_integrate_rmass_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void final_integrate_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void final_integrate_rmass_item(int) const;

 private:
  class AtomKokkos *atomKK;

  typename ArrayTypes<DeviceType>::t_x_array x;
  typename ArrayTypes<DeviceType>::t_v_array v;
  typename ArrayTypes<DeviceType>::t_f_array_const f;
  double *rmass;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread mass;
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

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/

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
FixStyle(deform/kk,FixDeformKokkos<LMPDeviceType>);
FixStyle(deform/kk/device,FixDeformKokkos<LMPDeviceType>);
FixStyle(deform/kk/host,FixDeformKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_DEFORM_KOKKOS_H
#define LMP_FIX_DEFORM_KOKKOS_H

#include "fix_deform.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType>
class FixDeformKokkos : public FixDeform {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixDeformKokkos(class LAMMPS *, int, char **);
  ~FixDeformKokkos();

  void init() override;
  void pre_exchange() override;
  void end_of_step() override;

 private:
  class DomainKokkos *domainKK;

  typename AT::t_x_array d_x;
  typename AT::t_int_1d d_mask;

  Kokkos::DualView<Set*, Kokkos::LayoutRight, LMPDeviceType> k_set;
  Kokkos::View<Set*,Kokkos::LayoutRight,DeviceType> d_set;

  KOKKOS_INLINE_FUNCTION
  void virtual apply_volume();

  KOKKOS_INLINE_FUNCTION
  void apply_strain();

  KOKKOS_INLINE_FUNCTION
  void update_domain();
};

}

#endif
#endif


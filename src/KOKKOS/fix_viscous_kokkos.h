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
FixStyle(viscous/kk,FixViscousKokkos<LMPDeviceType>);
FixStyle(viscous/kk/device,FixViscousKokkos<LMPDeviceType>);
FixStyle(viscous/kk/host,FixViscousKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_VISCOUS_KOKKOS_H
#define LMP_FIX_VISCOUS_KOKKOS_H

#include "fix_viscous.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixViscous{};

template<class DeviceType>
class FixViscousKokkos : public FixViscous {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixViscousKokkos(class LAMMPS *, int, char **);
  ~FixViscousKokkos() override;
  void init() override;
  void post_force(int) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixViscous, const int&) const;

 private:
  typename AT::t_v_array v;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d_randomread type;

  Kokkos::DualView<double*, Kokkos::LayoutRight, DeviceType> k_gamma;
};

}

#endif
#endif

